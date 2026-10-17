#!/usr/bin/env python3
"""XVI-R R03: G family (c*.G(c,h)) in the reduced system: is there a global coefficient vector on its pieces
(X_K: c* -> D^T Ktilde;  X_C: its cbar partner rho T;  X_S: h -> s0 h, antisymmetrised) that makes every G image
satisfy the reduced rows + ghost equation?  Exact nullspace over all 33 generators jointly."""
import json, pickle, sys  # noqa: E401
from pathlib import Path
import sympy as sp
R = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import r03_fix as F  # noqa: E402
S = F.S
OUT = R / "build" / "logs"


def main():
    d = pickle.loads((S.LOG / "q08_triv.pkl").read_bytes())
    B = pickle.loads((S.LOG / "q11_blocks.pkl").read_bytes())
    Ab = F.DM(pickle.loads((OUT / "r03_reduced_hhc.pkl").read_bytes()))
    AAb = Ab.vstack(F.DM(B["ccK"])).vstack(F.DM(B["ccb"])).vstack(F.DM(B["bose"])).vstack(F.DM(d["Ag"]))
    n = S.NDIM
    K, C, Sx = [], [], []
    for g in S.family_G():
        a = g["anti"]
        K.append([a[i] if i < S.OFF["C"] else 0 for i in range(n)])
        C.append([a[i] if S.OFF["C"] <= i < S.OFF["S"] else 0 for i in range(n)])
        Sx.append(g["field"])
    cols = []
    for P in (K, C, Sx):
        M = (AAb * F.DM(sp.Matrix(P).T)).to_Matrix()
        cols.append(M.reshape(M.rows * M.cols, 1))
    ns = F.DM(sp.Matrix.hstack(*cols)).nullspace().to_Matrix()
    res = {"pieces": ["X_K", "X_C(partner)", "X_S"], "nullspace": [[str(x) for x in ns.row(i)] for i in range(ns.rows)]}
    (OUT / "r03_gscan.json").write_text(json.dumps(res, indent=1) + "\n")
    print(res)


if __name__ == "__main__":
    main()
