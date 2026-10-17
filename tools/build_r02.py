#!/usr/bin/env python3
"""XVI-R R02 build: full R01 build (stage0 + R01 science + R01 gate), then R02 science, certificates, tests, R02.pdf,
R02 gate.  No sync.  Prints BUILD_OK."""
import subprocess, sys  # noqa: E401
from pathlib import Path
R = Path(__file__).resolve().parents[1]


def run(*a):
    r = subprocess.run([sys.executable, *a], cwd=R, capture_output=True, text=True)
    print(f"[{a[0]}] rc={r.returncode}", r.stdout[-400:], r.stderr[-1500:], flush=True)
    if r.returncode or (a[0] == "tools/build_r.py" and "BUILD_OK" not in r.stdout):
        raise SystemExit(f"FAIL {a}")


def main():
    run("tools/build_r.py")
    run("code/r02/r02_sci.py")
    run("code/r02/r02_certs.py")
    run("-m", "pytest", "-q", "tests_r02")
    for _ in range(3):
        r = subprocess.run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error", "R02.tex"], cwd=R / "manuscript",
                           capture_output=True, text=True)
        if r.returncode:
            raise SystemExit("FAIL tex\n" + r.stdout[-2000:])
    run("code/r02/release_r02.py")
    print("BUILD_OK", flush=True)


if __name__ == "__main__":
    main()
