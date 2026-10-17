#!/usr/bin/env python3
"""XVI-R R08 front 2: cohomology count in a Gram-reduced scheme.
The R06/R07 count is in the formal H/T algebra (hat dots independent symbols).  A Gram-reduced scheme imposes the
strictly-4D Gram identities  det[v_i . v_j] = 0  (5 hat vectors) on the barred subspace; such combinations become
evanescent (vanish in 4D) and are added to the evanescent counterterm space.
(i) enumeration: a 5x5 Gram determinant has 10 vector slots; in every X sector and every row the T-containing
    monomials carry <= 8 hat slots, so Gram identities act ONLY on pure-H monomials; the only sector with 10 hat slots
    is X_H (a1^2 a2^2 a3^2, momentum degree 4), where rows = cols = {a1,a2,a3,p1,p2} is forced (a repeated vector in
    rows or cols kills the det): exactly ONE Gram polynomial G (checked: G in the 4D basis, nonzero).
(ii) G is Bose symmetric; its sigma-bar image (hhc channel A: a1 -> sym(p1, c)) vanishes identically (p1 repeated);
    it does not enter the ghost rows or the ghost equation.
(iii) extended classification (hat rule on all external probe lines, complete Y as R07): coordinates n+1 (G added),
    trivial directions = (EV z, lambda) with FD z = lambda g (4D remainder proportional to G)."""
import json, sys, time  # noqa: E401
from pathlib import Path
import sympy as sp
from sympy import QQ_I
from sympy.polys.matrices import DomainMatrix
R = Path(__file__).resolve().parents[2]
for p in ("r07/code", "r06/code", "r05/code", "r03/code"):
    sys.path.insert(0, str(R / p))
import r07_reclass as RC  # noqa: E402
C6, F5, S, G, F6 = RC.C6, RC.F5, RC.S, RC.G, RC.F6
DM, rk, hcat = RC.DM, RC.rk, RC.hcat
U = S.U
I = sp.I
OUT = R / "build/logs"


def gram():
    v = ["a1", "a2", "a3", "p1", "p2"]
    M = sp.Matrix(5, 5, lambda i, j: U.H(v[i], v[j]))
    return sp.expand(M.det())


def main():
    t0 = time.time()
    import pickle
    S.SIG = sp.sympify(json.loads((S.LOG / "q07_hhc_tree.json").read_text())["sigma_F"])
    import triv  # noqa: F401
    S.triv = sys.modules["triv"]
    res = {}
    Gp = gram()
    M4 = {s: C6.enum(p, k, "4d") for s, (p, k) in C6.SECT.items()}
    w, rem = S.triv.coords(Gp, M4["H"], "H")
    res["G_nonzero"] = Gp != 0
    res["G_n_terms"] = len(sp.Add.make_args(Gp))
    res["G_in_4D_basis"] = rem == 0
    bose = all(sp.expand(Gp - S.R.remap(Gp, mp)) == 0 for mp in (
        {"p1": S.R.P2, "p2": S.R.P1, "a1": "a2", "a2": "a1"}, {"p1": S.R.P3, "p2": S.R.P2, "a1": "a3", "a3": "a1"}))
    res["G_bose_symmetric"] = bose
    A_img = S.R.strip_t(-I * S.R.polarise(sp.expand(Gp), "a1", "p1", "gc1"))
    res["G_sigma_image_zero"] = sp.expand(A_img) == 0
    # slot census: hat slots available in T-containing monomials of each X sector (Gram needs 10)
    cen = {}
    for s, (p, k) in C6.SECT.items():
        slots = sum(p.values()) + k
        cen[s] = {"total_slots": slots, "max_hat_slots_evanescent": slots - 2, "gram_possible_4D": slots >= 10}
    res["slot_census"] = cen
    # extended classification under the full hat rule
    d = pickle.loads((S.LOG / "q08_triv.pkl").read_bytes())
    B = pickle.loads((S.LOG / "q11_blocks.pkl").read_bytes())
    Ag = DM(d["Ag"])
    n = Ag.shape[1]
    off4, o = {}, 0
    for s in "KCSH":
        off4[s] = o
        o += len(M4[s])
    n4 = o
    gvec = [0] * n4
    for i, c in enumerate(w):
        gvec[off4["H"] + i] = c
    Ab = DM(pickle.loads((OUT / "r07_reduced_hhc_hat.pkl").read_bytes()))
    bk, bb = DM(RC.under("hatgen", RC.lin_block, "ccK")), DM(RC.under("hatgen", RC.lin_block, "ccb"))
    Abar = Ab.vstack(bk).vstack(bb).vstack(DM(B["bose"]))
    EV, FD, tags = RC.under("hatgen", RC.images, M4, n, n4, off4)
    EVm, FDm = DM(sp.Matrix(EV).T), DM(sp.Matrix(FD).T)
    # base (R07) classification
    Z = FDm.nullspace().transpose()
    Tev = EVm * Z
    cols = [list(Tev.to_Matrix()[:, j]) for j in range(Tev.shape[1])]
    AAb = Abar.vstack(Ag)
    cl0, Nb, T, Bg = F5.classify(AAb, Abar, Ag, cols, n)
    res["R07_reproduced"] = cl0
    # extended: coordinate n (G); its column in every row block is zero (image zero, not in ghost rows / ghost eq)
    zc = DomainMatrix.zeros((Abar.shape[0], 1), QQ_I)
    Abar_e = hcat(Abar, zc)
    Ag_e = hcat(Ag, DomainMatrix.zeros((Ag.shape[0], 1), QQ_I))
    AAb_e = Abar_e.vstack(Ag_e)
    FDg = hcat(FDm, DM(sp.Matrix([[-x] for x in gvec])))
    Zg = FDg.nullspace().transpose()
    EVg = EVm.vstack(DomainMatrix.zeros((1, EVm.shape[1]), QQ_I))
    EVg = hcat(EVg, DM(sp.Matrix([[0] * n + [1]]).T))
    Tg = EVg * Zg
    cols_e = [list(Tg.to_Matrix()[:, j]) for j in range(Tg.shape[1])]
    cl1, Nb1, T1, Bg1 = F5.classify(AAb_e, Abar_e, Ag_e, cols_e, n + 1)
    gdir = DM(sp.Matrix([[0] * n + [1]]).T)
    res["G_trivial"] = rk(hcat(T1, gdir)) == rk(T1)
    res["G_blind"] = rk(hcat(Bg1, gdir)) == rk(Bg1)
    res["gram_scheme"] = cl1
    res["H_formal"] = cl0["dim_H"]
    res["H_gram_reduced"] = cl1["dim_H"]
    res["seconds"] = round(time.time() - t0, 1)
    (OUT / "r08_gram.json").write_text(json.dumps(res, indent=1, default=str) + "\n")
    print(json.dumps(res, indent=1, default=str))


if __name__ == "__main__":
    main()
