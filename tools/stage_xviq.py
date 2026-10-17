#!/usr/bin/env python3
"""XVI-R stage 0: verify the pinned XVI-Q Q11 source-only zip, extract it into build/xviq and run XVI-Q's own Q11 build
there (no sync).  XVI-Q itself is never touched; its build products (certificates, build/logs) are the R inputs."""
import hashlib, shutil, subprocess, sys, zipfile  # noqa: E401
from pathlib import Path
R = Path(__file__).resolve().parents[1]
ZIP = R / "imports/XVI-Q-Q11-source-only.zip"
SHA = "1888f3277b793ee8cbf3954b3819bec03e41d8c4dc59757744dfda8b147a35ad"


def main():
    h = hashlib.sha256(ZIP.read_bytes()).hexdigest()
    if h != SHA:
        raise SystemExit(f"pinned XVI-Q zip mismatch {h}")
    dst = R / "build" / "xviq"
    if dst.exists():
        shutil.rmtree(dst)
    dst.mkdir(parents=True)
    zipfile.ZipFile(ZIP).extractall(dst)
    q = dst / "XVI-Q"
    r = subprocess.run([sys.executable, "tools/build_q06_q07_q08_q09_q10_q11.py"], cwd=q, capture_output=True, text=True)
    (R / "build" / "logs" / "stage_xviq.log").write_text(r.stdout + r.stderr)
    if r.returncode or "BUILD_OK" not in r.stdout:
        raise SystemExit("XVI-Q Q11 build from pinned zip failed")
    print("STAGE_XVIQ_OK")


if __name__ == "__main__":
    main()
