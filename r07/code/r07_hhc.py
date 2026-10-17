#!/usr/bin/env python3
"""XVI-R R07 front 1: the h.h.c loop residual in the reduced (Gamma-bar) variables.
R06 found: the stored XVI-Q hhc residual (channels A + B + C + D(EH+gf) + sigma_F E, full on-shell-B formulation)
paired with the REDUCED counterterm channels (A, D_EH) gives an unsolvable inhomogeneous fit.  In Gamma-bar, with
K~ = K + cbar F'(h), the gauge-fixing part of S_2 and the cbar channel E are absorbed into K~: the reduced row is
    Rbar = A + B2 + B3 + C + D2_EH + D3_EH        (loop channels; same graphs, gf piece of S_2 and E removed)
Checks: (t) reduced TREE row A_t + B_t(EH) (+ the E-less calibration) ; (l) loop identity full - reduced = D_gf + sigma E ;
(f) joint inhomogeneous fit: Rbar + reduced CT channels + hat-aware ccK/ccb rows + Bose + ghost equation.
Frozen code untouched (channel values read from frozen logs; D_EH recomputed with frozen functions)."""
import hashlib, json, sys, time  # noqa: E401
from pathlib import Path
import sympy as sp
R = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(R / "r06/code"))
import r06_refit as F6  # noqa: E402
B5, E2, G, S = F6.B5, F6.E2, F6.G, F6.S
hhc = S.hhc
I = sp.I
OUT = R / "build/logs"
LOG = S.LOG


def rd(k):
    return json.loads((LOG / f"q07_hhc_{k}.json").read_text())


def chan_D_split(i):
    ai, ki, _ = hhc.EXT[i]
    Gk = hhc._relabel(hhc.load_any("G03") + hhc.load_any("G07"), i)
    Gk = S.R.polarise(sp.expand(Gk), "aK", "KU", "KW")
    eh, gf = G.S2_EH_gf(hhc.ext(ai), G.probe_mat("E"), ki)
    return (hhc.strip(S.R.contract_id(Gk, I * eh, "K", "E", "h")), hhc.strip(S.R.contract_id(Gk, I * gf, "K", "E", "h")))


def tree_reduced():
    from hhc import TT, p1, p2, p3, GC, EXT  # noqa: F401
    Es = TT.sym(p1, GC, -I)
    A = G.V3_grav((Es, p1), (hhc.ext("a2"), p2), (hhc.ext("a3"), p3))
    Beh, Bgf, E = 0, 0, 0
    for i, j in (("2", "3"), ("3", "2")):
        ai, ki, _ = EXT[i]
        aj, kj, _ = EXT[j]
        M = G.Lc_h(hhc.ext(aj), kj, GC, p1)
        M.items = [(sp.Rational(1, 2) * c, u, w) for c, u, w in M.items]
        M.eta = sp.Rational(1, 2) * M.eta
        eh, gf = G.S2_EH_gf(hhc.ext(ai), M, ki)
        Beh += I * eh
        Bgf += I * gf
        E += G.V_cb_c_h(G.Dvec(hhc.ext(ai), ki), GC, p1, hhc.ext(aj), kj)
    A, Beh, Bgf, E = (hhc.strip(x) for x in (A, Beh, Bgf, E))
    sig = sp.sympify(rd("tree")["sigma_F"])
    return {"A_plus_B_EH_zero": sp.expand(A + Beh) == 0, "B_gf_plus_sigma_E_zero": sp.expand(Bgf + sig * E) == 0,
            "full_tree_zero": sp.expand(A + Beh + Bgf + sig * E) == 0}


def main():
    t0 = time.time()
    res = {"tree": tree_reduced()}
    sig = sp.sympify(rd("tree")["sigma_F"])
    P = {"A": sp.sympify(rd("A")["value"]), "C": sp.sympify(rd("C")["value"])}
    for i in ("2", "3"):
        b = rd("B" + i)
        for s in ("grav", "ghost", "scalar"):
            P[f"B{i}_{s}"] = sp.sympify(b[s])
    Dgf, Esum = 0, 0
    for i in ("2", "3"):
        deh, dgf = chan_D_split(i)
        P["D" + i + "_EH"] = deh
        Dgf += dgf
        Esum += sig * sp.sympify(rd("E" + i)["value"])
        res[f"D{i}_split_sums_to_frozen"] = sp.expand(deh + dgf - sp.sympify(rd("D" + i)["value"])) == 0
    Rbar = sp.expand(sum(P.values()))
    full = sp.sympify(json.loads((LOG / "q07_hhc_row.json").read_text())["residual"])
    diff = sp.expand(full - Rbar)
    gfE = sp.expand(Dgf + Esum)
    res["full_minus_reduced_equals_Dgf_plus_sigmaE"] = sp.expand(diff - gfE) == 0
    res["Dgf_plus_sigmaE_4D"] = str(F6.t0proj(gfE))
    res["Dgf_plus_sigmaE_4D_zero"] = F6.t0proj(gfE) == 0
    res["Dgf_plus_sigmaE_n_terms"] = len(sp.Add.make_args(gfE)) if gfE != 0 else 0
    res["Rbar_n_terms"] = len(sp.Add.make_args(Rbar)) if Rbar != 0 else 0
    res["Rbar_4D_zero"] = F6.t0proj(Rbar) == 0
    res["Rbar_sha256"] = hashlib.sha256(str(Rbar).encode()).hexdigest()
    (OUT / "r07_Rbar.txt").write_text(str(Rbar) + "\n")
    # joint fit, reduced + hat
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
    G.contract = B5.contract_hat
    try:
        rK = sp.expand(E2.ghost_row_total("ccK", XK, XS))
        rB = sp.expand(E2.ghost_row_total("ccb", XC, XS))
    finally:
        G.contract = B5._orig
    beq = []
    for mp in ({"p1": S.R.P2, "p2": S.R.P1, "a1": "a2", "a2": "a1"}, {"p1": S.R.P3, "p2": S.R.P2, "a1": "a3", "a3": "a1"}):
        beq += F6.coeffs(XH - S.R.remap(XH, mp), X)
    geqs = F6.coeffs(XC - S.RHO * S.Tmap(XK), X)
    eqs = F6.coeffs(rK, X) + F6.coeffs(rB, X) + F6.coeffs(sp.expand(Rbar + red), X) + beq
    A, b = sp.linear_eq_to_matrix(eqs, xall)
    Ag, bg = sp.linear_eq_to_matrix(geqs, xall)
    M, mb = A.col_join(Ag), b.col_join(bg)
    rM, rMb = F6.rk(M), F6.rk(M.row_join(mb))
    res["fit"] = {"n_eq": len(eqs) + len(geqs), "rank": rM, "rank_aug": rMb, "solvable_with_ghost_eq": rM == rMb,
                  "dim_N_with_ghost_eq": len(xall) - rM, "solvable_without_ghost_eq": F6.rk(A) == F6.rk(A.row_join(b))}
    if rM == rMb:
        sol = list(sp.linsolve(list(M * sp.Matrix(xall) - mb), xall))[0]
        free = set().union(*[sp.sympify(v).free_symbols for v in sol]) & X
        val = {x: sp.sympify(v).xreplace({f: 0 for f in free}) for x, v in zip(xall, sol)}
        cts = {n: sp.expand(e.xreplace(val)) for n, e in (("K", XK), ("C", XC), ("S", XS), ("H", XH))}
        res["fit"]["particular_CT_n_terms"] = {n: len(sp.Add.make_args(e)) if e != 0 else 0 for n, e in cts.items()}
        res["fit"]["particular_CT_all_evanescent"] = all(F6.t0proj(e) == 0 for e in cts.values())
        res["fit"]["verified_rows_zero"] = all(sp.expand(r.xreplace(val)) == 0 for r in (rK, rB, sp.expand(Rbar + red)))
        (OUT / "r07_Xbar.json").write_text(json.dumps({n: str(e) for n, e in cts.items()}, indent=1) + "\n")
    res["seconds"] = round(time.time() - t0, 1)
    (OUT / "r07_hhc.json").write_text(json.dumps(res, indent=1) + "\n")
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
