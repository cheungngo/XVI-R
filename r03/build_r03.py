#!/usr/bin/env python3
"""XVI-R R03 build: full R02 build (stage0, R01, R02 incl. gates), then the R03 chain, certificates, tests, R03.pdf, gate."""
import subprocess, sys  # noqa: E401
from pathlib import Path
R = Path(__file__).resolve().parents[1]


def run(*a):
    r = subprocess.run([sys.executable, *a], cwd=R, capture_output=True, text=True)
    print(f"[{a[0]}] rc={r.returncode}", r.stdout[-300:], r.stderr[-1500:], flush=True)
    if r.returncode or (a[0].startswith("tools/build") and "BUILD_OK" not in r.stdout):
        raise SystemExit(f"FAIL {a}")


def main():
    run("tools/build_r02.py")
    for s in ("r03_diag", "r03_legs", "r03_fix", "r03_fix2", "r03_final", "r03_gscan", "r03_certs"):
        run(f"r03/code/{s}.py")
    run("-m", "pytest", "-q", "r03/tests")
    for _ in range(3):
        r = subprocess.run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error", "R03.tex"], cwd=R / "r03",
                           capture_output=True, text=True)
        if r.returncode:
            raise SystemExit("FAIL tex\n" + r.stdout[-2000:])
    run("r03/code/release_r03.py")
    print("BUILD_OK", flush=True)


if __name__ == "__main__":
    main()
