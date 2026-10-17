#!/usr/bin/env python3
"""XVI-R R07: full re-run of the reduced classification and of the reduced inhomogeneous fit with the hat-aware
contraction on ALL external probe lines (c* line X<->ZZ as in R05, ghost line X<->Z, graviton line K~<->E as found in
R07, module r07_hatgen).  Two passes:
  'R05rule' : only the c*-line hat rule (reproduces R05/R06: sanity),
  'hatgen'  : hat rule on every external probe line.
Y = complete local space (evanescent + 4D generators, as R06), trivial evanescent directions = EV . null(FD).
Front 1 (reduced hhc loop residual Rbar = A + B + C + D_EH) re-solved under 'hatgen'."""
import json, pickle, sys, time  # noqa: E401
from pathlib import Path
import sympy as sp
from sympy import QQ_I
from sympy.polys.matrices import DomainMatrix
R = Path(__file__).resolve().parents[2]
for p in ("r07/code", "r06/code", "r05/code", "r03/code"):
    sys.path.insert(0, str(R / p))
import r07_hatgen as HG  # noqa: E402
import r06_complete as C6  # noqa: E402
import r07_hhc as H7  # noqa: E402
F6, B5, G, S, E2 = HG.F6, HG.B5, HG.G, HG.S, HG.F6.E2
F5 = C6.F5
DM, rk, hcat = C6.DM, C6.rk, C6.hcat
I = sp.I
OUT = R / "build/logs"
HG.EXTP |= {"Z"}
RULES = {"R05rule": B5.contract_hat, "hatgen": HG.contract_hatgen}


def under(rule, fn, *a):
    G.contract = RULES[rule]
    try:
        return fn(*a)
    finally:
        G.contract = B5._orig


def lin_block(target):
    XK, xk = F6.ansatz({"aK": 2, "gc1": 1, "a3": 2}, 3, "xk")
    XC, xc = F6.ansatz({"gb2": 1, "gc1": 1, "a3": 2}, 4, "xc")
    XS, xs = F6.ansatz({"z": 1, "gc1": 1, "gc2": 1}, 3, "xs")
    XH, xh = F6.ansatz({"a1": 2, "a2": 2, "a3": 2}, 4, "xh")
    xall = xk + xc + xs + xh
    XA = XK if target == "ccK" else XC
    r = sp.expand(E2.ghost_row_total(target, XA, XS) - E2.ghost_row_total(target, sp.Integer(0), sp.Integer(0)))
    gens = sorted([s for s in r.free_symbols if s not in set(xall)], key=str)
    return sp.linear_eq_to_matrix(sp.Poly(r, *gens).coeffs(), xall)[0]


def images(M4, n, n4, off4):
    EV, FD, tags = [], [], []
    for fam, (p, k) in C6.FAMS.items():
        for m in C6.enum(p, k, "all"):
            im = C6.images_4d(fam, m)
            if im is None:
                continue
            ev, fd = [0] * n, [0] * n4
            for sec, ex in im.items():
                if ex == 0:
                    continue
                v, w, rem = C6.split(ex, sec, M4)
                assert rem == 0
                for i, c in enumerate(v):
                    ev[S.OFF[sec] + i] += c
                for i, c in enumerate(w):
                    fd[off4[sec] + i] += c
            if any(sp.expand(x) != 0 for x in ev + fd):
                EV.append(ev)
                FD.append(fd)
                tags.append(fam)
    return EV, FD, tags


def main():
    t0 = time.time()
    S.SIG = sp.sympify(json.loads((S.LOG / "q07_hhc_tree.json").read_text())["sigma_F"])
    import triv  # noqa: F401
    S.triv = sys.modules["triv"]
    d = pickle.loads((S.LOG / "q08_triv.pkl").read_bytes())
    B = pickle.loads((S.LOG / "q11_blocks.pkl").read_bytes())
    Ag = DM(d["Ag"])
    n = Ag.shape[1]
    M4 = {s: C6.enum(p, k, "4d") for s, (p, k) in C6.SECT.items()}
    off4, o = {}, 0
    for s in "KCSH":
        off4[s] = o
        o += len(M4[s])
    lab = [(sec, mm) for sec in "KCSH" for mm in S.MS[sec]]
    idx = {(sec, F5.F4.tex_mono(mm)): i for i, (sec, mm) in enumerate(lab)}
    reps = []
    for b in json.loads((R / "certificates/R05_H_basis.json").read_text())["basis"]:
        v = [0] * n
        for sec, mm, c in b["terms"]:
            v[idx[(sec, mm)]] = F5.F4.coef(c)
        reps.append(v)
    C18 = DM(sp.Matrix(reps).T)
    res = {}
    store = {}
    for rule in ("R05rule", "hatgen"):
        t1 = time.time()
        if rule == "R05rule":
            Ab = DM(pickle.loads((OUT / "r03_reduced_hhc.pkl").read_bytes()))
            NB = pickle.loads((OUT / "r05_blocks.pkl").read_bytes())
            bk, bb = DM(NB["ccK"]), DM(NB["ccb"])
        else:
            Ab = DM(pickle.loads((OUT / "r07_reduced_hhc_hat.pkl").read_bytes()))
            bk, bb = DM(under(rule, lin_block, "ccK")), DM(under(rule, lin_block, "ccb"))
            res["ghost_blocks_equal_R05"] = {t: F5.same_rowspace(x, DM(NB0[t])) for t, x in (("ccK", bk), ("ccb", bb))}
        if rule == "R05rule":
            NB0 = NB
        Abar = Ab.vstack(bk).vstack(bb).vstack(DM(B["bose"]))
        AAb = Abar.vstack(Ag)
        EV, FD, tags = under(rule, images, M4, n, sum(len(v) for v in M4.values()), off4)
        EVm, FDm = DM(sp.Matrix(EV).T), DM(sp.Matrix(FD).T)
        Z = FDm.nullspace().transpose()
        Tev = EVm * Z
        cols = [list(Tev.to_Matrix()[:, j]) for j in range(Tev.shape[1])]
        cl, Nb, T, Bg = F5.classify(AAb, Abar, Ag, cols, n)
        W = hcat(T, Bg)
        cl["n_generators"] = {f: tags.count(f) for f in C6.FAMS}
        cl["dim_pure_evanescent_combos"] = Z.shape[1]
        cl["images_inside_Nbar"] = (AAb * Tev).is_zero_matrix
        cl["R05_18_in_Nbar"] = (AAb * C18).is_zero_matrix
        cl["R05_18_trivial_count"] = rk(C18) + rk(W) - rk(hcat(C18, W))
        cl["seconds"] = round(time.time() - t1, 1)
        res[rule] = cl
        store[rule] = (AAb, Abar, Nb, T, Bg)
        print(rule, cl, flush=True)
    # canonical basis of the hatgen cohomology (rref of Nb modulo W)
    AAb, Abar, Nb, T, Bg = store["hatgen"]
    W = hcat(T, Bg)
    basis = []
    cur = W
    for j in range(Nb.shape[1]):
        c = Nb.extract(list(range(n)), [j])
        if rk(hcat(cur, c)) > rk(cur):
            cur = hcat(cur, c)
            basis.append([str(x) for x in c.to_Matrix()])
    res["hatgen_H_basis_vectors"] = basis
    # front 1 under hatgen: Rbar + reduced CT channels + full ghost rows + Bose + ghost eq
    res["front1"] = under("hatgen", front1)
    res["seconds"] = round(time.time() - t0, 1)
    (OUT / "r07_reclass.json").write_text(json.dumps(res, indent=1, default=str) + "\n")
    print(json.dumps({k: v for k, v in res.items() if k != "hatgen_H_basis_vectors"}, indent=1, default=str))


def front1():
    hhc = S.hhc
    rd = H7.rd
    sig = sp.sympify(rd("tree")["sigma_F"])
    P = {"A": sp.sympify(rd("A")["value"]), "C": hhc.chan_C()}
    for i in ("2", "3"):
        b = rd("B" + i)
        for s in ("grav", "ghost", "scalar"):
            P[f"B{i}_{s}"] = sp.sympify(b[s])
        P["D" + i] = H7.chan_D_split(i)[0]
    Rbar = sp.expand(sum(P.values()))
    (OUT / "r07_Rbar_hatgen.txt").write_text(str(Rbar) + "\n")
    XK, xk = F6.ansatz({"aK": 2, "gc1": 1, "a3": 2}, 3, "xk")
    XC, xc = F6.ansatz({"gb2": 1, "gc1": 1, "a3": 2}, 4, "xc")
    XS, xs = F6.ansatz({"z": 1, "gc1": 1, "gc2": 1}, 3, "xs")
    XH, xh = F6.ansatz({"a1": 2, "a2": 2, "a3": 2}, 4, "xh")
    xall = xk + xc + xs + xh
    X = set(xall)
    red = S.R.strip_t(-I * S.R.polarise(sp.expand(XH), "a1", "p1", "gc1"))
    for i in ("2", "3"):
        ai, ki, Fi = hhc.EXT[i]
        Gk = S.R.polarise(sp.expand(hhc._relabel(XK, i)), "aK", "KU", "KW")
        eh, gf = G.S2_EH_gf(hhc.ext(ai), G.probe_mat("E"), ki)
        red += S.R.strip_t(S.R.contract_id(Gk, I * eh, "K", "E", "h"))
    rK = sp.expand(E2.ghost_row_total("ccK", XK, XS))
    rB = sp.expand(E2.ghost_row_total("ccb", XC, XS))
    out = {"Rbar_4D_zero": F6.t0proj(Rbar) == 0, "Rbar_n_terms": len(sp.Add.make_args(Rbar)),
           "ghost_rows_4D_zero": all(F6.t0proj(r.xreplace({x: 0 for x in xall})) == 0 for r in (rK, rB)),
           "sigma_F": str(sig)}
    beq = []
    for mp in ({"p1": S.R.P2, "p2": S.R.P1, "a1": "a2", "a2": "a1"}, {"p1": S.R.P3, "p2": S.R.P2, "a1": "a3", "a3": "a1"}):
        beq += F6.coeffs(XH - S.R.remap(XH, mp), X)
    geqs = F6.coeffs(XC - S.RHO * S.Tmap(XK), X)
    eqs = F6.coeffs(rK, X) + F6.coeffs(rB, X) + F6.coeffs(sp.expand(Rbar + red), X) + beq
    A, b = sp.linear_eq_to_matrix(eqs, xall)
    Ag, bg = sp.linear_eq_to_matrix(geqs, xall)
    M, mb = A.col_join(Ag), b.col_join(bg)
    rM, rMb = F6.rk(M), F6.rk(M.row_join(mb))
    out.update({"rank": rM, "rank_aug": rMb, "solvable_with_ghost_eq": rM == rMb, "dim_N_with_ghost_eq": len(xall) - rM,
                "solvable_without_ghost_eq": F6.rk(A) == F6.rk(A.row_join(b))})
    if rM == rMb:
        sol = list(sp.linsolve(list(M * sp.Matrix(xall) - mb), xall))[0]
        free = set().union(*[sp.sympify(v).free_symbols for v in sol]) & X
        val = {x: sp.sympify(v).xreplace({f: 0 for f in free}) for x, v in zip(xall, sol)}
        cts = {nm: sp.expand(e.xreplace(val)) for nm, e in (("K", XK), ("C", XC), ("S", XS), ("H", XH))}
        out["particular_CT_n_terms"] = {nm: len(sp.Add.make_args(e)) if e != 0 else 0 for nm, e in cts.items()}
        out["particular_CT_all_evanescent"] = all(F6.t0proj(e) == 0 for e in cts.values())
        out["verified_rows_zero"] = all(sp.expand(r.xreplace(val)) == 0 for r in (rK, rB, sp.expand(Rbar + red)))
        (OUT / "r07_Xbar_hatgen.json").write_text(json.dumps({nm: str(e) for nm, e in cts.items()}, indent=1) + "\n")
    return out


if __name__ == "__main__":
    main()
