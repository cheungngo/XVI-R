#!/usr/bin/env python3
"""XVI-R R03 diagnosis of the hhc sigma-composition gap.  The F image is split into channel pieces
  D_EH (K -> dS_EH/dh), D_gf (K -> dS_gf/dh), A_K (h -> s0 h on K.F), E (cbar-partner -> sigma_F F), A_C (h -> s0 h on the
partner), and G/P likewise; then we ask for ALL global coefficient vectors c with  A . sum_i c_i piece_i(Y) = 0  for every
generator Y simultaneously (a sign/normalisation/double-counting fix must be one of these).  Exact over QQ<I>."""
import json, pickle, sys, time  # noqa: E401
from pathlib import Path
import sympy as sp
from sympy import QQ_I
from sympy.polys.matrices import DomainMatrix
R = Path(__file__).resolve().parents[2]
QQ = R / "build" / "xviq" / "XVI-Q"
sys.path.insert(0, str(QQ / "code/p6q06/q07/q08/q09/q10/q11"))
import sigma11 as S  # noqa: E402
OUT = R / "build" / "logs"
I = sp.I


def D_parts(f):
    eh_t, gf_t = 0, 0
    for i in "123":
        fi = S.R.remap(f, S.REL[i]) if S.REL[i] else f
        fi = S.R.polarise(fi, "aK", "KU", "KW")
        eh, gf = S.G.S2_EH_gf(S.hhc.ext(f"a{i}"), S.G.probe_mat("E"), S.MOM[i])
        eh_t += S.R.strip_t(S.R.contract_id(fi, I * eh, "K", "E", "h"))
        gf_t += S.R.strip_t(S.R.contract_id(fi, I * gf, "K", "E", "h"))
    return eh_t, gf_t


def DM(M):
    return DomainMatrix.from_Matrix(M).convert_to(QQ_I)


def main():
    t0 = time.time()
    S.SIG = sp.sympify(json.loads((S.LOG / "q07_hhc_tree.json").read_text())["sigma_F"])
    d = pickle.loads((S.LOG / "q08_triv.pkl").read_bytes())
    A, Ag = DM(d["A"]), DM(d["Ag"])
    z = [0] * S.NDIM
    P = {k: [] for k in ("D_EH", "D_gf", "A_K", "E", "A_C")}
    check = []
    for m in S.monomials({"aK": 2, "a2": 2, "a3": 2}, 2):
        f = sp.expand(m + S.R.remap(m, S.SW23))
        if f == 0:
            continue
        phi = S.cbar_partner(f)
        eh, gf = D_parts(f)
        P["D_EH"].append(S.vec("H", eh)[0])
        P["D_gf"].append(S.vec("H", gf)[0])
        P["A_K"].append(S.vec("K", S.A_on_K(f))[0])
        P["E"].append(S.vec("H", S.E_op(phi))[0] if phi != 0 else z)
        P["A_C"].append(S.vec("C", S.A_on_C(phi))[0] if phi != 0 else z)
        check.append(S.add(S.vec("H", S.D_op(f))[0], [-x for x in S.add(P["D_EH"][-1], P["D_gf"][-1])]))
    res = {"n_F": len(P["E"]), "D_split_exact": all(sp.expand(x) == 0 for v in check for x in v)}
    names = list(P)
    blocks = {}
    for k in names:
        M = DM(sp.Matrix(P[k]).T)                     # n x nF
        blocks[k] = (A * M, Ag * M)
    # coefficient system: columns = channels, rows = all (row, generator) entries flattened
    def flat(Dm):
        return Dm.to_Matrix().reshape(Dm.shape[0] * Dm.shape[1], 1)
    for lab, use_g in (("rows", False), ("rows_and_ghost_eq", True)):
        cols = []
        for k in names:
            v = flat(blocks[k][0])
            if use_g:
                v = v.col_join(flat(blocks[k][1]))
            cols.append(v)
        C = DM(sp.Matrix.hstack(*cols))
        ns = C.nullspace().to_Matrix()
        res[f"channel_coefficients_{lab}"] = {"rank": C.rank(), "nullspace": [[str(x) for x in ns.row(i)] for i in range(ns.rows)]}
    # per-channel defect ranks and pairwise span relations
    res["rank_A_piece"] = {k: blocks[k][0].rank() for k in names}
    der = blocks["D_EH"][0] + blocks["D_gf"][0] + blocks["E"][0] - blocks["A_K"][0] - blocks["A_C"][0]
    res["rank_defect_derived"] = der.rank()
    noE = blocks["D_EH"][0] - blocks["A_K"][0] - blocks["A_C"][0]
    res["rank_defect_noGF_noE"] = noE.rank()
    res["rank_defect_EH_plusE"] = (blocks["D_EH"][0] + blocks["E"][0] - blocks["A_K"][0] - blocks["A_C"][0]).rank()
    res["rank_defect_full_noE"] = (blocks["D_EH"][0] + blocks["D_gf"][0] - blocks["A_K"][0] - blocks["A_C"][0]).rank()
    res["seconds"] = round(time.time() - t0, 1)
    (OUT / "r03_diag.json").write_text(json.dumps(res, indent=1) + "\n")
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
