#!/usr/bin/env python3
import subprocess, sys
from pathlib import Path
R = Path(__file__).resolve().parents[1]


def run(*a):
    r = subprocess.run([sys.executable, *a], cwd=R, capture_output=True, text=True)
    print(f"[{a[0]}] rc={r.returncode}", r.stdout[-300:], r.stderr[-1500:], flush=True)
    if r.returncode or (a[0].endswith("build_r07.py") and "BUILD_OK" not in r.stdout):
        raise SystemExit(f"FAIL {a}")


def main():
    run("r07/build_r07.py")
    for s in ("r08_gram", "r08_inventory", "r08_certs"):
        run(f"r08/code/{s}.py")
    run("-m", "pytest", "-q", "r08/tests")
    for _ in range(3):
        r = subprocess.run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error", "R08.tex"],
                           cwd=R / "r08", capture_output=True, text=True)
        if r.returncode:
            raise SystemExit("FAIL tex\n" + r.stdout[-2000:])
    run("r08/code/release_r08.py")
    print("BUILD_OK", flush=True)


if __name__ == "__main__":
    main()
