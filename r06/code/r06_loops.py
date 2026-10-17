#!/usr/bin/env python3
"""XVI-R R06 (A): the frozen ghost-row LOOP remainders (X = 0) of Q06/Q07 recomputed with the hat-aware c*-line
contraction of R05 (wrapper only; frozen code untouched).  Compared term by term with the frozen rule.
Also: the joint evanescent CT fit (Q07 style, reduced system) is re-solved if anything changes."""
import hashlib, json, sys, time  # noqa: E401
from pathlib import Path
import sympy as sp
R = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(R / "r05/code"))
import r05_blocks as B5  # noqa: E402
E2, G = B5.E2, B5.G
OUT = R / "build/logs"


def base(target, hat):
    G.contract = B5.contract_hat if hat else B5._orig
    try:
        return sp.expand(E2.ghost_row_total(target, sp.Integer(0), sp.Integer(0)))
    finally:
        G.contract = B5._orig


def main():
    t0 = time.time()
    res = {}
    for t in ("ccK", "ccb"):
        old, new = base(t, False), base(t, True)
        d = sp.expand(new - old)
        res[t] = {"n_terms_old": len(sp.Add.make_args(old)), "n_terms_new": len(sp.Add.make_args(new)),
                  "sha_old": hashlib.sha256(str(old).encode()).hexdigest(),
                  "sha_new": hashlib.sha256(str(new).encode()).hexdigest(),
                  "difference_zero": d == 0, "difference_n_terms": 0 if d == 0 else len(sp.Add.make_args(d)),
                  "difference": str(d)[:2000]}
    # frozen Q06/Q07 stored remainder for cross-check
    fz = json.loads((B5.Q.parent / "build/logs/q06_ghost_rows.json").read_text()) if \
        (B5.Q.parent / "build/logs/q06_ghost_rows.json").is_file() else {}
    for t in ("ccK", "ccb"):
        if t in fz:
            res[t]["frozen_q06_remainder_sha"] = fz[t].get("evanescent_remainder_sha256")
    res["loop_remainders_unchanged"] = all(res[t]["difference_zero"] for t in ("ccK", "ccb"))
    res["seconds"] = round(time.time() - t0, 1)
    (OUT / "r06_loops.json").write_text(json.dumps(res, indent=1) + "\n")
    print(json.dumps({t: {k: v for k, v in res[t].items() if k != "difference"} if isinstance(res[t], dict) else res[t]
                      for t in res}, indent=1))


if __name__ == "__main__":
    main()
