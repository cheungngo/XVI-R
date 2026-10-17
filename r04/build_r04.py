#!/usr/bin/env python3
import subprocess, sys
from pathlib import Path
R = Path(__file__).resolve().parents[1]


def run(*a):
    r = subprocess.run([sys.executable, *a], cwd=R, capture_output=True, text=True)
    print(f"[{a[0]}] rc={r.returncode}", r.stdout[-300:], r.stderr[-1500:], flush=True)
    if r.returncode or (a[0].endswith("build_r03.py") and "BUILD_OK" not in r.stdout):
        raise SystemExit(f"FAIL {a}")


def main():
    run("r03/build_r03.py")
    for s in ("r04_gderive", "r04_gpat", "r04_s0hat", "r04_final", "r04_certs"):
        run(f"r04/code/{s}.py")
    run("-m", "pytest", "-q", "r04/tests")
    for _ in range(3):
        r = subprocess.run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error", "R04.tex"],
                           cwd=R / "r04", capture_output=True, text=True)
        if r.returncode:
            raise SystemExit("FAIL tex\n" + r.stdout[-2000:])
    run("r04/code/release_r04.py")
    print("BUILD_OK", flush=True)


if __name__ == "__main__":
    main()
