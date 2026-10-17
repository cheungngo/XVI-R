#!/usr/bin/env python3
"""XVI-R R06 (A) propagation: re-solve the joint evanescent CT fit with the hat-aware contraction applied to the WHOLE
ccK/ccb ghost rows (loop remainder + linear part); frozen code untouched (contract swapped by wrapper only).
Fits:  (Q) Q07/Q08 system: full hhc row (stored residual, channels A,D,E) + ccK + ccb + Bose (+ ghost equation);
       (R) reduced Gamma-bar system: stored hhc residual + reduced channels (A, D_EH) + ccK + ccb + Bose (+ ghost eq).
Also: 4D projection (T -> 0) of new vs old remainders, i.e. the Q07 '4D closure' statement."""
import json, sys, time  # noqa: E401
from pathlib import Path
import sympy as sp
from sympy.polys.matrices import DomainMatrix
R = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(R / "r05/code"))
sys.path.insert(0, str(R / "r03/code"))
import r05_blocks as B5  # noqa: E402
E2, G = B5.E2, B5.G
import r03_fix as F3  # noqa: E402
S = F3.S
ansatz = B5.ansatz
OUT = R / "build/logs"
I = sp.I


def t0proj(e):
    return sp.expand(sp.sympify(e).xreplace({s: 0 for s in sp.sympify(e).free_symbols if s.name.startswith("T_")}))


def coeffs(r, X):
    r = sp.expand(r)
    if r == 0:
        return []
    gens = sorted([s for s in r.free_symbols if s not in X], key=str)
    return sp.Poly(r, *gens).coeffs()


def rk(M):
    return DomainMatrix.from_Matrix(M).convert_to(sp.QQ_I).rank() if M.shape[0] else 0


def main():
    t0 = time.time()
    res = {}
    XK, xk = ansatz({"aK": 2, "gc1": 1, "a3": 2}, 3, "xk")
    XC, xc = ansatz({"gb2": 1, "gc1": 1, "a3": 2}, 4, "xc")
    XS, xs = ansatz({"z": 1, "gc1": 1, "gc2": 1}, 3, "xs")
    XH, xh = ansatz({"a1": 2, "a2": 2, "a3": 2}, 4, "xh")
    xall = xk + xc + xs + xh
    X = set(xall)
    sig = sp.sympify(json.loads((S.LOG / "q07_hhc_tree.json").read_text())["sigma_F"])
    rows = {}
    for hat in (False, True):
        G.contract = B5.contract_hat if hat else B5._orig
        try:
            rows[hat] = {"ccK": sp.expand(E2.ghost_row_total("ccK", XK, XS)),
                         "ccb": sp.expand(E2.ghost_row_total("ccb", XC, XS))}
            base = {t: sp.expand(E2.ghost_row_total(t, sp.Integer(0), sp.Integer(0))) for t in ("ccK", "ccb")}
        finally:
            G.contract = B5._orig
        for t in ("ccK", "ccb"):
            res.setdefault("loop_4D_projection", {})[f"{t}_{'hat' if hat else 'old'}_T0"] = str(t0proj(base[t]))
    for t in ("ccK", "ccb"):
        res["loop_4D_projection"][f"{t}_4D_part_unchanged"] = \
            res["loop_4D_projection"][f"{t}_hat_T0"] == res["loop_4D_projection"][f"{t}_old_T0"]
    hhc_full = E2.hhc_ansatz_total(XH, XK, XC, sig)
    base_h = sp.sympify(json.loads((S.LOG / "q07_hhc_row.json").read_text())["residual"])
    red = S.R.strip_t(-I * S.R.polarise(sp.expand(XH), "a1", "p1", "gc1"))
    for i in ("2", "3"):
        ai, ki, Fi = S.hhc.EXT[i]
        Gk = S.R.polarise(sp.expand(S.hhc._relabel(XK, i)), "aK", "KU", "KW")
        eh, gf = S.G.S2_EH_gf(S.hhc.ext(ai), S.G.probe_mat("E"), ki)
        red += S.R.strip_t(S.R.contract_id(Gk, I * eh, "K", "E", "h"))
    hhc_red = sp.expand(base_h + red)
    beq = []
    for mp in ({"p1": S.R.P2, "p2": S.R.P1, "a1": "a2", "a2": "a1"}, {"p1": S.R.P3, "p2": S.R.P2, "a1": "a3", "a3": "a1"}):
        beq += coeffs(XH - S.R.remap(XH, mp), X)
    geqs = coeffs(XC - S.RHO * S.Tmap(XK), X)
    for sysname, hrow in (("Q_full_hhc", hhc_full), ("R_reduced_hhc", hhc_red)):
        for hat in (False, True):
            eqs = coeffs(rows[hat]["ccK"], X) + coeffs(rows[hat]["ccb"], X) + coeffs(hrow, X) + beq
            A, b = sp.linear_eq_to_matrix(eqs, xall)
            Ag, bg = sp.linear_eq_to_matrix(geqs, xall)
            rA, rAb = rk(A), rk(A.row_join(b))
            M, mb = A.col_join(Ag), b.col_join(bg)
            rM, rMb = rk(M), rk(M.row_join(mb))
            key = f"{sysname}_{'hat' if hat else 'old'}"
            res[key] = {"n_eq": len(eqs), "rank": rA, "solvable": rA == rAb, "dim_N": len(xall) - rA,
                        "rank_with_ghost_eq": rM, "solvable_with_ghost_eq": rM == rMb, "dim_N_with_ghost_eq": len(xall) - rM}
            if hat and rM == rMb:
                sol = list(sp.linsolve(list(M * sp.Matrix(xall) - mb), xall))[0]
                free = set().union(*[sp.sympify(v).free_symbols for v in sol]) & X
                val = {x: sp.sympify(v).xreplace({f: 0 for f in free}) for x, v in zip(xall, sol)}
                cts = {n: sp.expand(e.xreplace(val)) for n, e in (("K", XK), ("C", XC), ("S", XS), ("H", XH))}
                res[key]["particular_CT_n_terms"] = {n: len(sp.Add.make_args(e)) if e != 0 else 0 for n, e in cts.items()}
                res[key]["particular_CT_all_evanescent"] = all(t0proj(e) == 0 for e in cts.values())
                res[key]["verified_rows_zero"] = all(sp.expand(r.xreplace(val)) == 0 for r in
                                                     (rows[True]["ccK"], rows[True]["ccb"], hrow))
                res[key]["particular_CT"] = {n: str(e)[:3000] for n, e in cts.items()}
    res["seconds"] = round(time.time() - t0, 1)
    (OUT / "r06_refit.json").write_text(json.dumps(res, indent=1) + "\n")
    print(json.dumps({k: ({kk: vv for kk, vv in v.items() if kk not in ("particular_CT",) and not kk.endswith("T0")}
                          if isinstance(v, dict) else v) for k, v in res.items()}, indent=1))


if __name__ == "__main__":
    main()
