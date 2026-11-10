"""Chunked multipart upload to a Zenodo deposition's files bucket.

Why this exists: a single-stream PUT of a 500 MB file over a long-lived TCP
connection is fragile — any 502, nginx hiccup or client-side timeout kills
the whole upload and we restart from byte 0. Multipart lets us upload in
small (20 MB) independently-retriable chunks. If chunk 7 fails, we retry
just chunk 7.

Zenodo supports multipart via the Invenio-Files-REST API:
  1. POST {bucket}/{key}?uploads=&size=TOTAL&part_size=CHUNK
  2. For each chunk N (1..last): PUT {bucket}/{key}?uploadId=X&partNumber=N
  3. POST {bucket}/{key}?uploadId=X  (to complete / commit)
Reference: https://inveniordm.docs.cern.ch/reference/rest_api_drafts/
Zenodo also supports "normal" PUT, which this script also supports as a
fallback for files smaller than CHUNK_SIZE.

Usage:
  export ZENODO_TOKEN=...
  python3 zenodo_upload_chunked.py <deposit_id> <file> [<file> ...]
"""
from __future__ import annotations
import argparse
import hashlib
import os
import sys
import time
from pathlib import Path

import requests

ZENODO_API = "https://zenodo.org/api"
CHUNK_SIZE = 20 * 1024 * 1024   # 20 MB per chunk
RETRIES_PER_CHUNK = 10
CONNECT_TIMEOUT = 30


def log(m):
    print(f"[chunked] {m}", flush=True)


def _auth(token):
    return {"Authorization": f"Bearer {token}"}


def _retry(fn, retries=RETRIES_PER_CHUNK, label=""):
    for i in range(1, retries + 1):
        try:
            r = fn()
            if 200 <= r.status_code < 300:
                return r
            log(f"  {label} attempt {i}/{retries}: HTTP {r.status_code}  body={r.text[:200]}")
        except (requests.exceptions.ConnectionError,
                requests.exceptions.ChunkedEncodingError,
                requests.exceptions.ReadTimeout,
                requests.exceptions.Timeout) as e:
            log(f"  {label} attempt {i}/{retries}: {type(e).__name__}: {str(e)[:200]}")
        if i < retries:
            backoff = min(10 + 10 * i, 180)
            log(f"  wait {backoff}s then retry")
            time.sleep(backoff)
    raise RuntimeError(f"{label} failed after {retries} attempts")


def get_bucket(deposit_id, token):
    r = _retry(
        lambda: requests.get(f"{ZENODO_API}/deposit/depositions/{deposit_id}",
                             headers=_auth(token), timeout=CONNECT_TIMEOUT),
        label="GET deposit",
    )
    return r.json()["links"]["bucket"]


def list_existing(deposit_id, token):
    r = _retry(
        lambda: requests.get(f"{ZENODO_API}/deposit/depositions/{deposit_id}/files",
                             headers=_auth(token), timeout=CONNECT_TIMEOUT),
        label="GET files",
    )
    return {f["filename"]: f for f in r.json()}


def delete_existing(deposit_id, token, file_id):
    _retry(
        lambda: requests.delete(
            f"{ZENODO_API}/deposit/depositions/{deposit_id}/files/{file_id}",
            headers=_auth(token), timeout=CONNECT_TIMEOUT),
        label=f"DELETE file {file_id}",
    )


def md5_of(path, chunk=1024 * 1024):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def upload_simple(bucket, path, token):
    """Fallback for files < CHUNK_SIZE: single PUT, Content-Type set correctly."""
    size = path.stat().st_size
    log(f"  simple PUT ({size/1024/1024:.1f} MB)")
    headers = {**_auth(token),
               "Content-Type": "application/octet-stream",
               "Content-Length": str(size)}
    def do():
        with open(path, "rb") as f:
            return requests.put(f"{bucket.rstrip('/')}/{path.name}",
                                data=f, headers=headers, timeout=(CONNECT_TIMEOUT, None))
    r = _retry(do, label=f"PUT {path.name}")
    return r.json()


def upload_multipart(bucket, path, token):
    """Chunked multipart: initiate -> N PUTs -> complete."""
    size = path.stat().st_size
    n_parts = (size + CHUNK_SIZE - 1) // CHUNK_SIZE
    log(f"  multipart: {size/1024/1024:.1f} MB in {n_parts} parts of up to {CHUNK_SIZE/1024/1024:.0f} MB")

    # 1. Initiate
    init_url = f"{bucket.rstrip('/')}/{path.name}?uploads=&size={size}&partSize={CHUNK_SIZE}"
    r = _retry(
        lambda: requests.post(init_url, headers={**_auth(token), "Content-Type": "application/json"},
                              timeout=CONNECT_TIMEOUT),
        label="initiate multipart",
    )
    j = r.json()

    # 2. Upload parts
    t0 = time.time()
    with open(path, "rb") as f:
        for n in range(1, n_parts + 1):
            buf = f.read(CHUNK_SIZE)
            actual_len = len(buf)
            def do(buf=buf, actual_len=actual_len, n=n):
                return requests.put(
                    f"{bucket.rstrip('/')}/{path.name}?uploadId={j.get('key', path.name)}&partNumber={n}",
                    data=buf,
                    headers={**_auth(token),
                             "Content-Type": "application/octet-stream",
                             "Content-Length": str(actual_len)},
                    timeout=(CONNECT_TIMEOUT, 600),
                )
            _retry(do, label=f"part {n}/{n_parts}")
            elapsed = time.time() - t0
            uploaded_mb = n * CHUNK_SIZE / 1024 / 1024
            if uploaded_mb > size / 1024 / 1024:
                uploaded_mb = size / 1024 / 1024
            rate = uploaded_mb / max(elapsed, 0.01)
            log(f"  part {n}/{n_parts} done  ({uploaded_mb:.1f} MB sent, avg {rate:.2f} MB/s)")

    # 3. Complete
    r = _retry(
        lambda: requests.post(f"{bucket.rstrip('/')}/{path.name}?uploadId={j.get('key', path.name)}",
                              headers=_auth(token), timeout=CONNECT_TIMEOUT),
        label="complete multipart",
    )
    return r.json()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("deposit_id")
    ap.add_argument("files", nargs="+", type=Path)
    ap.add_argument("--no-overwrite", action="store_true")
    args = ap.parse_args()

    token = os.environ.get("ZENODO_TOKEN")
    if not token:
        sys.exit("set ZENODO_TOKEN")

    log(f"deposit {args.deposit_id}")
    bucket = get_bucket(args.deposit_id, token)
    log(f"bucket: {bucket}")
    existing = list_existing(args.deposit_id, token)
    log(f"already uploaded: {len(existing)} files")

    for p in args.files:
        log(f"=== {p.name} ({p.stat().st_size/1024/1024:.1f} MB) ===")
        if not p.exists():
            log(f"  SKIP: file not found")
            continue
        if p.name in existing:
            if args.no_overwrite:
                log(f"  SKIP: already present (--no-overwrite)")
                continue
            log(f"  deleting existing (id={existing[p.name]['id']})")
            delete_existing(args.deposit_id, token, existing[p.name]["id"])

        local_md5 = md5_of(p)
        log(f"  local md5: {local_md5}")

        if p.stat().st_size <= CHUNK_SIZE:
            result = upload_simple(bucket, p, token)
        else:
            try:
                result = upload_multipart(bucket, p, token)
            except Exception as e:
                log(f"  multipart failed ({e}), falling back to simple PUT")
                result = upload_simple(bucket, p, token)

        remote_md5 = (result.get("checksum") or "").removeprefix("md5:")
        if remote_md5 and remote_md5 == local_md5:
            log(f"  OK checksums match")
        elif remote_md5:
            log(f"  WARN: remote md5 {remote_md5} != local {local_md5}")
        else:
            log(f"  (no checksum returned in response; verify via file listing)")

    log("done.")


if __name__ == "__main__":
    sys.exit(main() or 0)
