#!/usr/bin/env python3
"""XVI-R build: stage0 (XVI-Q Q11 from pinned zip) -> R01 science -> certificates -> tests -> note PDF -> release gate.
No sync anywhere.  Prints BUILD_OK."""
import shutil, subprocess, sys  # noqa: E401
from pathlib import Path
R = Path(__file__).resolve().parents[1]


def run(*a):
    r = subprocess.run([sys.executable, *a], cwd=R, capture_output=True, text=True)
    print(f"[{a[0]}] rc={r.returncode}", r.stdout[-400:], r.stderr[-1500:], flush=True)
    if r.returncode:
        raise SystemExit(f"FAIL {a}")


def main():
    for d in ("build", "certificates"):
        shutil.rmtree(R / d, ignore_errors=True)
    (R / "build/logs").mkdir(parents=True)
    run("tools/stage_xviq.py")
    run("code/r01/r01_sci.py")
    run("code/r01/r01_certs.py")
    run("-m", "pytest", "-q", "tests")
    for _ in range(2):
        r = subprocess.run(["pdflatex", "-interaction=nonstopmode", "R01.tex"], cwd=R / "manuscript", capture_output=True, text=True)
    if r.returncode:
        raise SystemExit("FAIL tex\n" + r.stdout[-2000:])
    run("code/r01/release_r01.py")
    print("BUILD_OK", flush=True)


if __name__ == "__main__":
    main()
