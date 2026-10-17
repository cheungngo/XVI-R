#!/usr/bin/env python3
"""XVI-R R01 science on top of the imported XVI-Q Q11 (built from the pinned zip in build/xviq).
Goal 1 (characterise the hhc sigma-composition gap): for every K~-based F / G sigma-image that is NOT in N, ask in
which coordinate sector a correction delta can restore A(v+delta)=0, with delta restricted to one sector at a time
(K, C=cbar, S=c*cc, H=hhh) and to hat degree 1 / >=2 coordinates (every coordinate is evanescent); report solvability counts, the rank
of the defect span, and whether the corrected vector also satisfies the ghost equation.  A solvable correction is a
LOCALISATION of the defect, not a derivation of a sigma-image (no direction is thereby certified trivial).
Goal 2 (the 43 unclassified directions of N_g): invariant 'visibility' q(S) = dim N_g / (W + N_g cap {x_S = 0}),
W = sigma-span(passing) + blind-in-N_g, for coordinate sets S (sectors; hat degree).  q(S)=0 means every
unclassified direction has a representative vanishing on S.  Classification itself (trivial/physical) is NOT claimed."""
import json, pickle, sys, time  # noqa: E401
from pathlib import Path
import sympy as sp
from sympy.polys.matrices import DomainMatrix
from sympy import QQ_I
R = Path(__file__).resolve().parents[2]
QQ = R / "build" / "xviq" / "XVI-Q"
sys.path.insert(0, str(QQ / "code/p6q06/q07/q08/q09/q10/q11"))
import sigma11 as S  # noqa: E402
OUT = R / "build" / "logs"


def rk(M):
    return DomainMatrix.from_Matrix(M).rank() if M.shape[1] and M.shape[0] else 0


def coords():
    lab = {}
    for sec in "KCSH":
        for i, m in enumerate(S.MS[sec]):
            tdeg = sum(e for b, e in sp.sympify(m).as_powers_dict().items() if getattr(b, "name", "").startswith("T_"))
            lab[S.OFF[sec] + i] = (sec, tdeg)
    return lab


def DM(M):
    return DomainMatrix.from_Matrix(M).convert_to(QQ_I)


def left_kernel(M):
    """DomainMatrix whose rows span {L : L*M = 0} (exact, QQ<I>)."""
    if M.shape[1] == 0:
        return DM(sp.eye(M.shape[0]))
    return DM(M.T).nullspace()


def in_colspace(L, b):
    return L.shape[0] == 0 or (L * b).is_zero_matrix


def main():
    t0 = time.time()
    S.SIG = sp.sympify(json.loads((S.LOG / "q07_hhc_tree.json").read_text())["sigma_F"])
    d = pickle.loads((S.LOG / "q08_triv.pkl").read_bytes())
    A, Ag = d["A"], d["Ag"]
    n = A.shape[1]
    lab = coords()
    sets = {f"sector_{s}": [j for j in range(n) if lab[j][0] == s] for s in "KCSH"}
    # every ansatz coordinate is evanescent (>= 1 hatted T factor); split by hat degree instead
    sets["all_coords"] = list(range(n))
    sets["hat_deg_1"] = [j for j in range(n) if lab[j][1] == 1]
    sets["hat_deg_ge2"] = [j for j in range(n) if lab[j][1] >= 2]
    sets["CS"] = sets["sector_C"] + sets["sector_S"]
    res = {"n_coords": n, "set_sizes": {k: len(v) for k, v in sets.items()}}
    g1 = {}
    LK = {}
    AA = A.col_join(Ag)
    for name, fam in (("F", S.family_F()), ("G", S.family_G())):
        fails = []
        for g in fam:
            v = S.add(g["anti"], S.scale(-1, g["field"]))
            if any(sp.expand(x) != 0 for x in v) and any(sp.expand(x) != 0 for x in A * sp.Matrix(v)):
                fails.append(v)
        st = {"n_fail": len(fails)}
        if fails:
            Dm = sp.Matrix.hstack(*[A * sp.Matrix(v) for v in fails])
            st["rank_defect_span"] = rk(Dm)
            for k, idx in sets.items():
                if k not in LK:
                    LK[k] = (left_kernel(A[:, idx]), left_kernel(AA[:, idx]))
                LA, LAA = LK[k]
                ok = sum(in_colspace(LA, DM(A * sp.Matrix(v))) for v in fails)
                okg = sum(in_colspace(LAA, DM(AA * sp.Matrix(v))) for v in fails)
                st[k] = {"correctable": ok, "correctable_and_ghost_eq": okg,
                         "defect_span_inside_image": in_colspace(LA, DM(Dm))}
        g1[name] = st
    res["goal1_defect_localisation"] = g1
    keep = pickle.loads((S.LOG / "q11_sigma_vecs.pkl").read_bytes())
    V = keep[("F", "derived_minus1")] + keep[("G", "derived_minus1")]
    Mv = sp.Matrix.hstack(*[sp.Matrix(v) for v in V])
    full = A.col_join(Ag)
    NgM = DM(full).nullspace().transpose().to_Matrix()
    Ng = list(range(NgM.shape[1]))
    zc = [j for j in range(n) if all(A[i, j] == 0 for i in range(A.shape[0]))]
    Bm = sp.Matrix.hstack(*[sp.Matrix([1 if i == j else 0 for i in range(n)]) for j in zc])
    ker = [DM(Ag * Bm).nullspace().transpose().to_Matrix()[:, i] for i in range(DM(Ag * Bm).nullspace().shape[0])]
    W = Mv.row_join(sp.Matrix.hstack(*[Bm * k for k in ker]))
    rW = rk(W)
    vis = {}
    for k, idx in sets.items():
        P = sp.Matrix.vstack(*[NgM[j, :] for j in idx])
        kP = DM(P).nullspace().transpose().to_Matrix() if P.shape[0] else sp.eye(len(Ng))
        sub = NgM * kP if kP.shape[1] else sp.zeros(n, 0)
        vis[k] = len(Ng) - rk(W.row_join(sub))
    res["goal2"] = {"dim_N_g": len(Ng), "dim_W": rW, "unclassified": len(Ng) - rW, "visibility_q": vis}
    res["seconds"] = round(time.time() - t0, 1)
    (OUT / "r01_sci.json").write_text(json.dumps(res, indent=1, default=str) + "\n")
    print(json.dumps(res, indent=1, default=str))


if __name__ == "__main__":
    main()
