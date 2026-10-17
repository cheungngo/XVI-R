#!/usr/bin/env python3
"""Which G generators fail the reduced system (Q11 construction, eps=-1), and under which blocks? Pattern by monomial."""
import json, pickle, sys
from pathlib import Path
import sympy as sp
R = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(R / "r03/code"))
import r03_fix as F
S = F.S
OUT = R / "build" / "logs"


def main():
    S.SIG = sp.sympify(json.loads((S.LOG / "q07_hhc_tree.json").read_text())["sigma_F"])
    d = pickle.loads((S.LOG / "q08_triv.pkl").read_bytes())
    B = pickle.loads((S.LOG / "q11_blocks.pkl").read_bytes())
    Ab = F.DM(pickle.loads((OUT / "r03_reduced_hhc.pkl").read_bytes()))
    blk = {"hhc_red": Ab, "ccK": F.DM(B["ccK"]), "ccb": F.DM(B["ccb"]), "bose": F.DM(B["bose"]), "geq": F.DM(d["Ag"])}
    AAb = Ab.vstack(blk["ccK"]).vstack(blk["ccb"]).vstack(blk["bose"]).vstack(blk["geq"])
    TR = {"H_a3_a3", "T_a3_a3"}
    rows = []
    for m in S.monomials({"gv": 1, "gc1": 1, "a3": 2}, 2):
        # Q11 image
        xk = sp.expand(-sp.I * S.U.H("aK", "p2") * S.hhc.subst_vec(m, "gv", {"aK": 1}))
        xc = sp.expand(S.RHO * S.Tmap(sp.nsimplify(xk, rational=True))) if xk != 0 else 0
        s = S.R.remap(m, {"p1": S.R.P1, "p2": S.R.P3, "gv": "z"})
        s = S.R.strip_t(-sp.I * S.R.polarise(sp.expand(s), "a3", "p2", "gc2"))
        xs = sp.expand(s - S.R.remap(s, {"p1": S.R.P2, "p2": S.R.P1, "gc1": "gc2", "gc2": "gc1"}))
        v = S.add(S.vec("K", xk)[0], S.vec("C", xc)[0] if xc != 0 else [0]*S.NDIM, S.scale(-1, S.vec("S", xs)[0]))
        dv = F.DM(sp.Matrix(v))
        fail = [b for b, M in blk.items() if not (M * dv).is_zero_matrix]
        names = {s.name for s in m.free_symbols}
        rows.append({"m": str(m), "trace_h": bool(names & TR), "T_gc_gv": "T_gc1_gv" in names,
                     "fail": fail, "ok": not fail})
    res = {"n_ok": sum(r["ok"] for r in rows), "n_fail": sum(not r["ok"] for r in rows),
           "fail_trace": sum(not r["ok"] and r["trace_h"] for r in rows),
           "ok_trace": sum(r["ok"] and r["trace_h"] for r in rows),
           "fail_Tgc": sum(not r["ok"] and r["T_gc_gv"] for r in rows),
           "by_fail": {}}
    for r in rows:
        k = "PASS" if r["ok"] else " ".join(r["fail"])
        res["by_fail"][k] = res["by_fail"].get(k, 0) + 1
    res["failures"] = [r for r in rows if not r["ok"]]
    (OUT / "r04_gpat.json").write_text(json.dumps(res, indent=1) + "\n")
    print(json.dumps({k: res[k] for k in res if k != "failures"}, indent=1))
    for r in res["failures"]:
        print("FAIL", r["m"], r["fail"])


if __name__ == "__main__":
    main()
