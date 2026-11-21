"""Upload large files to a Zenodo draft deposition via the API.

Browser uploads for multi-GB files are flaky because browsers don't resume
partial uploads. The Zenodo REST API supports direct PUT uploads to the files
bucket, which is much more reliable.

Usage (interactive, run from a terminal):

    export ZENODO_TOKEN="your_token_here"   # generate at zenodo.org/account/settings/applications/tokens/new/
    python zenodo_upload.py <deposit_id> <file_path> [<file_path> ...]

The deposit ID is the number in your Zenodo draft URL
(zenodo.org/uploads/12345678 → 12345678).

Features:
- Deletes an existing file with the same name before re-uploading (clears stuck placeholders).
- Uses requests streaming so the whole file doesn't load into RAM.
- Prints a tqdm-style progress bar.
- Retries on connection errors with exponential backoff.

Example:
    python zenodo_upload.py 12345678 \\
        data/matrix/pan_malignant_matrix_v4.parquet \\
        data/matrix/enhancer_candidate_matrix.parquet
"""
from __future__ import annotations
import argparse
import os
import sys
import time
import hashlib
from pathlib import Path

import requests

ZENODO_API = "https://zenodo.org/api"


def log(msg: str) -> None:
    print(f"[zenodo-upload] {msg}", flush=True)


def _auth_headers(token: str) -> dict:
    """Pass the token as an Authorization header rather than a URL query parameter.
    Zenodo supports both but the header form is more reliable (not logged in URLs,
    not stripped by some proxies) and some endpoints 500 on the query-param form."""
    return {"Authorization": f"Bearer {token}"}


def _get_with_retry(url: str, token: str, retries: int = 4) -> requests.Response:
    for attempt in range(1, retries + 1):
        try:
            r = requests.get(url, headers=_auth_headers(token), timeout=30)
            if r.status_code < 500:
                return r
            log(f"  server {r.status_code} on attempt {attempt}/{retries}: {r.text[:200]}")
        except requests.exceptions.RequestException as e:
            log(f"  network error on attempt {attempt}/{retries}: {e}")
        if attempt < retries:
            backoff = 5 * (2 ** (attempt - 1))
            log(f"  retrying GET in {backoff}s")
            time.sleep(backoff)
    r.raise_for_status()
    return r


def get_bucket_url(deposit_id: str, token: str) -> str:
    r = _get_with_retry(f"{ZENODO_API}/deposit/depositions/{deposit_id}", token)
    r.raise_for_status()
    data = r.json()
    bucket = data.get("links", {}).get("bucket")
    if not bucket:
        raise RuntimeError(
            f"deposition {deposit_id} has no files bucket — is it already published? "
            f"API response keys: {list(data.get('links', {}).keys())}"
        )
    return bucket


def list_files(deposit_id: str, token: str) -> list[dict]:
    r = _get_with_retry(f"{ZENODO_API}/deposit/depositions/{deposit_id}/files", token)
    r.raise_for_status()
    return r.json()


def delete_file(deposit_id: str, file_id: str, token: str) -> None:
    r = requests.delete(f"{ZENODO_API}/deposit/depositions/{deposit_id}/files/{file_id}",
                        headers=_auth_headers(token), timeout=30)
    if r.status_code not in (204, 200):
        log(f"  delete returned {r.status_code}: {r.text[:200]}")
    r.raise_for_status()


class _TqdmReader:
    """Wrap a file-like object so requests can stream it while we show progress."""
    def __init__(self, path: Path, size: int):
        self.f = open(path, "rb")
        self.total = size
        self.sent = 0
        self.start = time.time()
        self._last_print = 0.0

    def read(self, chunk: int = -1) -> bytes:
        data = self.f.read(chunk if chunk and chunk > 0 else 1024 * 1024)
        self.sent += len(data)
        now = time.time()
        if now - self._last_print > 1.0 or len(data) == 0:
            self._print()
            self._last_print = now
        return data

    def _print(self) -> None:
        pct = 100 * self.sent / max(self.total, 1)
        elapsed = max(time.time() - self.start, 0.001)
        mb_s = (self.sent / 1024 / 1024) / elapsed
        bar = int(40 * self.sent / max(self.total, 1))
        sys.stdout.write(
            f"\r  [{'=' * bar}{' ' * (40 - bar)}] "
            f"{self.sent / 1024 / 1024:7.1f} / {self.total / 1024 / 1024:7.1f} MB "
            f"({pct:5.1f}%) at {mb_s:5.2f} MB/s"
        )
        sys.stdout.flush()

    def close(self):
        self.f.close()
        print()

    def __len__(self):
        return self.total


def put_file(bucket_url: str, path: Path, token: str, retries: int = 10) -> dict:
    size = path.stat().st_size
    attempt = 0
    while attempt < retries:
        attempt += 1
        log(f"  attempt {attempt}/{retries}: PUT {path.name} ({size / 1024 / 1024:.1f} MB)")
        reader = _TqdmReader(path, size)
        try:
            r = requests.put(
                f"{bucket_url.rstrip('/')}/{path.name}",
                data=reader,
                headers={**_auth_headers(token), "Content-Length": str(size)},
                timeout=(30, None),  # 30s connect, no read timeout (streaming)
            )
            reader.close()
            if 200 <= r.status_code < 300:
                return r.json()
            log(f"  server returned {r.status_code}: {r.text[:400]}")
        except (requests.exceptions.ConnectionError, requests.exceptions.ChunkedEncodingError, requests.exceptions.ReadTimeout) as e:
            reader.close()
            log(f"  network error: {e}")
        if attempt < retries:
            # Backoff tuned for the two failure modes we see:
            #   - client-side TCP timeout / connection reset → retry FAST (30-120s)
            #     because the stall was just a network blip, Zenodo itself is fine
            #   - Zenodo-side 502 nginx error → retry SLOWER (90s-5min) to let their
            #     server recover
            # Linear-ish growth with a 5-minute ceiling. Over 10 attempts this gives
            # a worst-case 10×5=50 min per file (vs 10×15=150 min with the old curve).
            backoff = min(30 + 30 * attempt, 300)
            log(f"  retrying in {backoff}s (= {backoff // 60}min {backoff % 60}s)")
            time.sleep(backoff)
    raise RuntimeError(f"failed to upload {path.name} after {retries} attempts")


def md5_file(path: Path, chunk: int = 1024 * 1024) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        while True:
            data = f.read(chunk)
            if not data: break
            h.update(data)
    return h.hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("deposit_id")
    ap.add_argument("files", nargs="+", type=Path)
    ap.add_argument("--no-overwrite", action="store_true",
                    help="if the file already exists in the deposit, skip (default: delete + re-upload)")
    args = ap.parse_args()

    token = os.environ.get("ZENODO_TOKEN")
    if not token:
        log("ERROR: set ZENODO_TOKEN env var with your personal access token")
        log("  export ZENODO_TOKEN='your_token_from_zenodo.org/account/settings/applications/tokens'")
        return 2

    for p in args.files:
        if not p.exists():
            log(f"ERROR: file not found: {p}")
            return 2

    log(f"deposit {args.deposit_id}: fetching bucket URL")
    bucket_url = get_bucket_url(args.deposit_id, token)
    log(f"  bucket: {bucket_url}")

    log("listing existing files in deposit")
    existing = {f["filename"]: f for f in list_files(args.deposit_id, token)}
    for name, meta in existing.items():
        log(f"  existing: {name}  ({meta.get('filesize', '?')} bytes, checksum={meta.get('checksum', '?')[:16]}…)")

    for p in args.files:
        log(f"=== {p.name} ===")
        if p.name in existing:
            if args.no_overwrite:
                log(f"  already exists, skipping (--no-overwrite)")
                continue
            log(f"  deleting existing (file_id={existing[p.name]['id']})")
            delete_file(args.deposit_id, existing[p.name]["id"], token)
        local_md5 = md5_file(p)
        log(f"  local md5: {local_md5}")
        result = put_file(bucket_url, p, token)
        remote_md5 = result.get("checksum", "").removeprefix("md5:")
        log(f"  remote md5: {remote_md5}")
        if local_md5 == remote_md5:
            log(f"  ✓ checksums match")
        else:
            log(f"  ✗ checksum MISMATCH — delete and retry")
            return 3
        log(f"  file size: {result.get('size')} bytes")

    log("all uploads complete")
    return 0


if __name__ == "__main__":
    sys.exit(main())
