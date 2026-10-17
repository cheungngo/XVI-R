#!/usr/bin/env python3
"""R04: apply s0 hat{h} = 0.  For Y=c*.G(c,h), the piece X_S = -s0h.dY/dh only receives the unhatted
projection of dY/dh.  Operationally: if G's dependence on h is purely hatted (every a3 index sits in a T_),
X_S := 0; more refined: polarise only the H_a3_* part of G (drop T_a3_* before polarise)."""
import json, pickle, sys
from pathlib import Path
import sympy as sp
R = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(R / "r03/code"))
import r03_fix as F
S = F.S
OUT = R / "build" / "logs"


def split_H_T(m, a="a3"):
    """m = mH + mT where mT has every term containing T_a_*."""
    m = sp.expand(m)
    mH, mT = 0, 0
    for t in sp.Add.make_args(m):
        if any(s.name.startswith("T_") and a in s.name.split("_") for s in t.free_symbols):
            mT += t
        else:
            mH += t
    return sp.expand(mH), sp.expand(mT)


def XS_from(m, mode):
    if mode == "full":
        src = m
    elif mode == "H_only":
        src, _ = split_H_T(m)
    elif mode == "zero_if_any_T_a3":
        _, mT = split_H_T(m)
        src = 0 if mT != 0 else m
    else:
        raise ValueError(mode)
    if src == 0:
        return 0
    s = S.R.remap(src, {"p1": S.R.P1, "p2": S.R.P3, "gv": "z"})
    s = S.R.strip_t(-sp.I * S.R.polarise(sp.expand(s), "a3", "p2", "gc2"))
    return sp.expand(s - S.R.remap(s, {"p1": S.R.P2, "p2": S.R.P1, "gc1": "gc2", "gc2": "gc1"}))


def main():
    S.SIG = sp.sympify(json.loads((S.LOG / "q07_hhc_tree.json").read_text())["sigma_F"])
    d = pickle.loads((S.LOG / "q08_triv.pkl").read_bytes())
    B = pickle.loads((S.LOG / "q11_blocks.pkl").read_bytes())
    Ab = F.DM(pickle.loads((OUT / "r03_reduced_hhc.pkl").read_bytes()))
    AAb = Ab.vstack(F.DM(B["ccK"])).vstack(F.DM(B["ccb"])).vstack(F.DM(B["bose"])).vstack(F.DM(d["Ag"]))
    res = {}
    for mode in ("full", "H_only", "zero_if_any_T_a3"):
        fails = 0
        n = 0
        for m in S.monomials({"gv": 1, "gc1": 1, "a3": 2}, 2):
            xk = sp.expand(-sp.I * S.U.H("aK", "p2") * S.hhc.subst_vec(m, "gv", {"aK": 1}))
            xc = sp.expand(S.RHO * S.Tmap(sp.nsimplify(xk, rational=True))) if xk != 0 else 0
            xs = XS_from(m, mode)
            v = S.add(S.vec("K", xk)[0], S.vec("C", xc)[0] if xc != 0 else [0] * S.NDIM,
                      S.scale(-1, S.vec("S", xs)[0] if xs != 0 else [0] * S.NDIM))
            if not any(sp.expand(x) != 0 for x in v):
                continue
            n += 1
            if not (AAb * F.DM(sp.Matrix(v))).is_zero_matrix:
                fails += 1
        res[mode] = {"n_images": n, "failing": fails}
    (OUT / "r04_s0hat.json").write_text(json.dumps(res, indent=1) + "\n")
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
