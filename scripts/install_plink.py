#!/usr/bin/env python
"""Locate or install PLINK 1.9 for this workspace.

Search order:
1. ``--dest`` directory (default: ``<workspace>/tools/plink``)
2. ``PLINK_EXE`` environment variable
3. ``plink.exe`` / ``plink`` on ``PATH``
4. download from the official S3 bucket

The official bucket host is DNS-poisoned on some networks.  When a direct
download fails, this script resolves the bucket through DNS-over-HTTPS and
downloads the archive with parallel HTTP Range requests (each chunk retried),
then verifies the zip and extracts the binary.

Usage:
    python install_plink.py [--dest DIR] [--version 20260927] [--threads 8]
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
import threading
import time
import urllib.request
import zipfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

HOST = "s3.amazonaws.com"
# 官方 S3 安装包按平台区分（2026-11-07 验证 20260927 三个包均可下载）。
ASSET_TEMPLATES = {
    "win32": "plink1-assets/plink_win64_{version}.zip",
    "linux": "plink1-assets/plink_linux_x86_64_{version}.zip",
    "darwin": "plink1-assets/plink_mac_{version}.zip",
}


def asset_path(version: str) -> str:
    """Return the official S3 asset path for this platform."""
    platform = sys.platform
    if platform.startswith("win"):
        key = "win32"
    elif platform.startswith("linux"):
        key = "linux"
    elif platform == "darwin":
        key = "darwin"
    else:
        raise RuntimeError(f"不支持的平台 {platform!r}：请手动下载 PLINK 1.9 后用 --dest 指定")
    return ASSET_TEMPLATES[key].format(version=version)
DOH_ENDPOINTS = [
    "https://doh.pub/dns-query?name={host}&type=A",
    "https://dns.alidns.com/resolve?name={host}&type=A",
]
FALLBACK_IPS = [
    "54.231.128.0", "16.182.69.0", "52.217.85.254", "16.15.165.220",
    "16.182.98.40", "16.15.212.180", "54.231.224.184", "16.15.237.37",
]
DEFAULT_VERSION = "20260927"


def log(msg: str) -> None:
    print(f"[install_plink] {msg}", flush=True)


def workspace_root() -> Path:
    # <workspace>/.agents/skills/plink-analysis/scripts/install_plink.py
    return Path(__file__).resolve().parents[4]


def find_existing(explicit_dest: Path | None) -> Path | None:
    candidates: list[Path] = []
    if explicit_dest:
        candidates.append(explicit_dest / ("plink.exe" if os.name == "nt" else "plink"))
    env = os.environ.get("PLINK_EXE")
    if env:
        candidates.append(Path(env))
    candidates.append(workspace_root() / "tools" / "plink" / ("plink.exe" if os.name == "nt" else "plink"))
    which = shutil.which("plink.exe" if os.name == "nt" else "plink")
    if which:
        candidates.append(Path(which))
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return None


def verify(exe: Path) -> str:
    try:
        proc = subprocess.run([str(exe), "--version"], capture_output=True, text=True,
                              timeout=60, encoding="utf-8", errors="replace")
        text = (proc.stdout or proc.stderr).strip().splitlines()
        return text[0] if text else "unknown"
    except Exception as exc:  # noqa: BLE001
        return f"无法运行: {exc}"


def resolve_ips() -> list[str]:
    ips: list[str] = []
    for tpl in DOH_ENDPOINTS:
        try:
            req = urllib.request.Request(tpl.format(host=HOST),
                                         headers={"accept": "application/dns-json"})
            with urllib.request.urlopen(req, timeout=15) as resp:
                data = json.load(resp)
            for answer in data.get("Answer", []):
                if answer.get("type") == 1:
                    ips.append(answer.get("data", ""))
        except Exception as exc:  # noqa: BLE001
            log(f"DoH 解析失败 ({tpl.split('/')[2]}): {exc}")
    seen: set[str] = set()
    out: list[str] = []
    for ip in ips + FALLBACK_IPS:
        if ip and ip.count(".") == 3 and ip not in seen:
            seen.add(ip)
            out.append(ip)
    return out


def _curl() -> str:
    return "curl.exe" if os.name == "nt" else "curl"


def content_length(url: str, ips: list[str]) -> int | None:
    for ip in ips:
        proc = subprocess.run([_curl(), "-sSI", "--max-time", "20",
                               "--resolve", f"{HOST}:443:{ip}", url],
                              capture_output=True, text=True)
        for line in proc.stdout.splitlines():
            if line.lower().startswith("content-length:"):
                return int(line.split(":", 1)[1].strip())
    return None


def download(url: str, total: int, ips: list[str], threads: int, chunk: int,
             workdir: Path, retries: int = 8) -> Path:
    ranges = [(start, min(start + chunk - 1, total - 1)) for start in range(0, total, chunk)]
    log(f"下载 {total} 字节，{len(ranges)} 个分块，{threads} 线程")
    lock = threading.Lock()
    state = {"bytes": 0, "chunks": 0}

    def fetch(idx: int) -> None:
        start, end = ranges[idx]
        out = workdir / f"chunk_{idx:05d}.bin"
        if out.exists() and out.stat().st_size == end - start + 1:
            with lock:
                state["bytes"] += out.stat().st_size
                state["chunks"] += 1
            return
        last_error: Exception | None = None
        for attempt in range(retries):
            ip = ips[(idx + attempt) % len(ips)]
            tmp = workdir / f"chunk_{idx:05d}.tmp"
            tmp.unlink(missing_ok=True)
            proc = subprocess.run([_curl(), "-sS", "--fail", "--max-time", "120",
                                   "--resolve", f"{HOST}:443:{ip}",
                                   "-r", f"{start}-{end}", "-o", str(tmp), url],
                                  capture_output=True, text=True)
            if proc.returncode == 0 and tmp.exists() and tmp.stat().st_size == end - start + 1:
                tmp.replace(out)
                with lock:
                    state["bytes"] += out.stat().st_size
                    state["chunks"] += 1
                    if state["chunks"] % 10 == 0:
                        log(f"  {state['bytes']}/{total} 字节（{state['chunks']}/{len(ranges)} 块）")
                return
            last_error = RuntimeError(proc.stderr.strip()[:150] or "短读")
            tmp.unlink(missing_ok=True)
            time.sleep(0.5)
        raise RuntimeError(f"分块 {idx} 下载失败: {last_error}")

    with ThreadPoolExecutor(max_workers=threads) as pool:
        futures = [pool.submit(fetch, idx) for idx in range(len(ranges))]
        for future in as_completed(futures):
            future.result()

    archive = workdir / "plink.zip"
    with archive.open("wb") as handle:
        for idx in range(len(ranges)):
            handle.write((workdir / f"chunk_{idx:05d}.bin").read_bytes())
    return archive


def install(dest: Path, version: str, threads: int, chunk_kb: int) -> Path:
    dest.mkdir(parents=True, exist_ok=True)
    url = f"https://{HOST}/{asset_path(version)}"
    ips = resolve_ips()
    log(f"解析到 {len(ips)} 个 S3 IP：{', '.join(ips[:4])} ...")
    total = content_length(url, ips)
    if not total:
        raise RuntimeError("无法获取安装包大小（网络可能完全不可用）；请手动下载 PLINK 1.9 并解压到 "
                           f"{dest}，或用 --dest/PLINK_EXE 指定已有可执行文件")
    log(f"安装包大小: {total} 字节")
    workdir = Path(tempfile.mkdtemp(prefix="plink_install_"))
    try:
        archive = download(url, total, ips, threads, chunk_kb * 1024, workdir)
        with zipfile.ZipFile(archive) as zf:
            names = zf.namelist()
            member = next((n for n in names if n.lower().endswith(".exe") or n.lower() == "plink"), None)
            if member is None:
                raise RuntimeError(f"压缩包内未找到可执行文件: {names}")
            zf.extract(member, dest)
            extracted = dest / member
            target = dest / ("plink.exe" if os.name == "nt" else "plink")
            if extracted != target and extracted.is_file():
                shutil.move(str(extracted), str(target))
    finally:
        shutil.rmtree(workdir, ignore_errors=True)
    exe = dest / ("plink.exe" if os.name == "nt" else "plink")
    if os.name != "nt":
        # zipfile 提取不保留 Unix 可执行位，手动补上，否则无法运行
        os.chmod(exe, 0o755)
    log(f"已安装: {exe}")
    return exe


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dest", default=None, help="安装目录（默认 <workspace>/tools/plink）")
    parser.add_argument("--version", default=DEFAULT_VERSION)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--chunk-kb", type=int, default=256)
    parser.add_argument("--check-only", action="store_true")
    args = parser.parse_args()

    dest = Path(args.dest).expanduser().resolve() if args.dest else workspace_root() / "tools" / "plink"
    existing = find_existing(dest)
    if existing:
        log(f"已找到 PLINK: {existing}")
        version = verify(existing)
        log(f"版本: {version}")
        print(json.dumps({"ok": True, "exe": str(existing), "version": version}, ensure_ascii=False))
        return 0
    if args.check_only:
        log("未找到 PLINK")
        print(json.dumps({"ok": False, "exe": None, "version": None}, ensure_ascii=False))
        return 1
    exe = install(dest, args.version, args.threads, args.chunk_kb)
    version = verify(exe)
    log(f"版本: {version}")
    print(json.dumps({"ok": True, "exe": str(exe), "version": version}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
