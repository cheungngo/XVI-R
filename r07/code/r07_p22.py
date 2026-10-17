#!/usr/bin/env python3
"""XVI-R R07 front 2: p22 kappa^2 evanescent counterterms (XVI-P p22 EV-CT) versus the kappa^3 evanescent system.
(a) identification: the p22 EV-CTs are -(gt-structure part) of the kappa^2 two-point UV poles.  In XVI-Q's own
    calculus the gt-structure part of a pole o is  P_ev(o) = o - lift(o|_{T=0}),  lift: H(x,y) -> H(x,y)+T(x,y) for
    polarisation/probe pairs (g4 = ghat + gt; 4D structures are g4/p monomials).  Checked per graph against the pinned
    XVI-P certificates with the Q05 validation ratio (own/XVI-P) and the p22 coefficients (= -UV pole).
(b) cubic completion: insert X2 = -P_ev(two-point pole) for Gamma_hh, Gamma_K~c (D13), Gamma_cbar c (D06) into the
    kappa^3 rows (hhc channels B, C; ghost rows K, c2), all genuine three-point loops set to zero, hat rule on all
    external probe lines; solve  sigma-bar X3' + Xi = 0  with X3' in the evanescent ansatz (reduced system, Bose,
    ghost equation).  Solvable <=> the p22 operators admit an evanescent kappa^3 completion; X3' is defined modulo Nbar_g.
Frozen code untouched (module attributes swapped for the duration of a call)."""
import json, sys, time  # noqa: E401
from pathlib import Path
import sympy as sp
R = Path(__file__).resolve().parents[2]
for p in ("r07/code", "r06/code", "r05/code", "r03/code"):
    sys.path.insert(0, str(R / p))
import r07_hatgen as HG  # noqa: E402
F6, B5, G, S, E2 = HG.F6, HG.B5, HG.G, HG.S, HG.F6.E2
HG.EXTP |= {"Z"}
hhc = S.hhc
RW = S.R
I = sp.I
OUT = R / "build/logs"
QC = R / "build/xviq/XVI-Q/code"
sys.path.insert(0, str(QC / "p6q05"))
import validate as VA  # noqa: E402
PV = R / "build/xviq/XVI-Q/provenance/XVI-P/certificates"
MOMS = {"p1", "p2", "l"}


def lift(e):
    rep = {}
    for s in e.free_symbols:
        p = s.name.split("_")
        if len(p) == 3 and p[0] == "H" and p[1] not in MOMS and p[2] not in MOMS:
            rep[s] = s + sp.Symbol(f"T_{p[1]}_{p[2]}")
    return sp.expand(e.xreplace(rep))


def P_ev(o):
    o = sp.expand(o)
    return sp.expand(o - lift(F6.t0proj(o)))


def ratio(a, b):
    if b == 0:
        return {"ratio": None, "constant": a == 0}
    r = sp.cancel(a / b)
    return {"ratio": str(r), "constant": not r.free_symbols}


def ident():
    out = {}
    p22 = json.loads((PV / "P22_evanescent.json").read_text())["EV_CT"]
    refs = {"D01": json.loads((PV / "P16_batchC.json").read_text())["graphs"]["D01"]["result"],
            "D03": json.loads((PV / "P16_batchC.json").read_text())["graphs"]["D03"]["result"],
            "D11": json.loads((PV / "P18_D11.json").read_text())["result"]}
    g = hhc.gamma_hh(hhc.ext("a1"), hhc.p1, hhc.ext("a2"))
    own = {"D01": RW.strip_t(g[0]), "D03": RW.strip_t(g[1]), "D11": RW.strip_t(g[2])}
    for k in ("D01", "D03", "D11"):
        ref = refs[k]
        full = VA.hh_reference(ref)
        evref = VA.hh_reference({s: v for s, v in ref.items() if v.get("evanescent")})
        ct = p22["coefficients_hh_counterterm_per_graph"][k]
        ct_is_minus_pole = all(sp.simplify(sp.sympify(ct.get(s, 0)) + sp.sympify(v["UV_pole"])) == 0
                               for s, v in ref.items() if v.get("evanescent")) and \
            all(ref.get(s, {}).get("evanescent") for s in ct)
        rf, re_ = ratio(own[k], full), ratio(P_ev(own[k]), evref)
        out[k] = {"own_over_P_full": rf, "own_ev_over_P_ev": re_, "p22_CT_equals_minus_pole": ct_is_minus_pole,
                  "same_ratio": rf["constant"] and re_["constant"] and rf["ratio"] == re_["ratio"],
                  "own_ev_n_terms": len(sp.Add.make_args(P_ev(own[k])))}
    # D13 (K~.c, E2) and D06 (cbar c, SE) with the Q05 validation conventions
    p, l = U_p1, None
    d13 = sp.sympify(json.loads(json.dumps(VA.D13()))["mine_d4"])
    r13 = json.loads((PV / "P15_D13.json").read_text())["result"]
    H, T = VA.U.H, VA.U.T
    pp, aKp = H("p1", "p1"), H("aK", "p1")
    e2ref = sp.sympify(r13["E2"]["UV_pole"]) * 2 * T("aK", "gc") * aKp * pp
    out["D13"] = {"own_ev_over_P_E2": ratio(P_ev(d13), e2ref), "own_over_P_full": {k: v for k, v in VA.D13().items() if k != "mine_d4"},
                  "p22_CT_equals_minus_pole": sp.sympify(p22["coefficient_Kc_E2"]) == -sp.sympify(r13["E2"]["UV_pole"])}
    d06 = sp.sympify(VA.D06()["mine_d4"])
    D = json.loads((R / "build/xviq/XVI-Q/certificates/Q02_bulk_graphs.json").read_text())["D06"]["laurent"]
    seref = sp.sympify(D["SE_Gt^rho_nu_p4"]["UV_pole"]) * T("gb", "gc") * pp ** 2
    out["D06"] = {"own_ev_over_P_SE": ratio(P_ev(d06), seref), "own_over_P_full": {k: v for k, v in VA.D06().items() if k != "mine_d4"},
                  "p22_CT_equals_minus_pole": sp.sympify(p22["coefficient_ccb_SE"]) == -sp.sympify(D["SE_Gt^rho_nu_p4"]["UV_pole"])}
    return out


U_p1 = None


def completion():
    orig = {"gh": hhc.gamma_hh, "d13": RW.D13_own, "d13k": RW.D13_K_ext, "d06": RW.D06_own, "load": RW.load,
            "la": hhc.load_any}

    def X2hh(E1, k, E2_):
        g = orig["gh"](E1, k, E2_)
        return (-P_ev(sum(g)), sp.Integer(0), sp.Integer(0))
    try:
        hhc.gamma_hh = X2hh
        RW.D13_own = lambda *a: -P_ev(orig["d13"](*a))
        RW.D13_K_ext = lambda *a: -P_ev(RW.strip_t(orig["d13k"](*a)))
        RW.D06_own = lambda *a: -P_ev(orig["d06"](*a))
        RW.load = lambda g: sp.Integer(0)
        hhc.load_any = lambda g: sp.Integer(0)
        G.contract = HG.contract_hatgen
        Xi = {"B": sp.expand(sum(sum(hhc.chan_B(i).values()) for i in ("2", "3"))), "C": sp.expand(hhc.chan_C())}
        XK, xk = F6.ansatz({"aK": 2, "gc1": 1, "a3": 2}, 3, "xk")
        XC, xc = F6.ansatz({"gb2": 1, "gc1": 1, "a3": 2}, 4, "xc")
        XS, xs = F6.ansatz({"z": 1, "gc1": 1, "gc2": 1}, 3, "xs")
        XH, xh = F6.ansatz({"a1": 2, "a2": 2, "a3": 2}, 4, "xh")
        xall = xk + xc + xs + xh
        X = set(xall)
        rK = sp.expand(E2.ghost_row_total("ccK", XK, XS))
        rB = sp.expand(E2.ghost_row_total("ccb", XC, XS))
        red = S.R.strip_t(-I * S.R.polarise(sp.expand(XH), "a1", "p1", "gc1"))
        for i in ("2", "3"):
            ai, ki, Fi = hhc.EXT[i]
            Gk = S.R.polarise(sp.expand(hhc._relabel(XK, i)), "aK", "KU", "KW")
            eh, gf = G.S2_EH_gf(hhc.ext(ai), G.probe_mat("E"), ki)
            red += S.R.strip_t(S.R.contract_id(Gk, I * eh, "K", "E", "h"))
    finally:
        hhc.gamma_hh, RW.D13_own, RW.D13_K_ext, RW.D06_own, RW.load, hhc.load_any = (
            orig["gh"], orig["d13"], orig["d13k"], orig["d06"], orig["load"], orig["la"])
        G.contract = B5._orig
    z = {x: 0 for x in xall}
    xi_hhc = sp.expand(Xi["B"] + Xi["C"])
    out = {"Xi_terms": {"hhc": len(sp.Add.make_args(xi_hhc)) if xi_hhc != 0 else 0,
                        "ccK": len(sp.Add.make_args(rK.xreplace(z))) if rK.xreplace(z) != 0 else 0,
                        "ccb": len(sp.Add.make_args(rB.xreplace(z))) if rB.xreplace(z) != 0 else 0},
           "Xi_4D_zero": all(F6.t0proj(e) == 0 for e in (xi_hhc, rK.xreplace(z), rB.xreplace(z)))}
    beq = []
    for mp in ({"p1": S.R.P2, "p2": S.R.P1, "a1": "a2", "a2": "a1"}, {"p1": S.R.P3, "p2": S.R.P2, "a1": "a3", "a3": "a1"}):
        beq += F6.coeffs(XH - S.R.remap(XH, mp), X)
    geqs = F6.coeffs(XC - S.RHO * S.Tmap(XK), X)
    eqs = F6.coeffs(rK, X) + F6.coeffs(rB, X) + F6.coeffs(sp.expand(xi_hhc + red), X) + beq
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
        out["completion_n_terms"] = {nm: len(sp.Add.make_args(e)) if e != 0 else 0 for nm, e in cts.items()}
        out["completion_all_evanescent"] = all(F6.t0proj(e) == 0 for e in cts.values())
        out["verified_rows_zero"] = all(sp.expand(r.xreplace(val)) == 0 for r in (rK, rB, sp.expand(xi_hhc + red)))
        (OUT / "r07_p22_completion.json").write_text(json.dumps({nm: str(e) for nm, e in cts.items()}, indent=1) + "\n")
    return out


def main():
    t0 = time.time()
    S.SIG = sp.sympify(json.loads((S.LOG / "q07_hhc_tree.json").read_text())["sigma_F"])
    res = {"identification": ident()}
    print(json.dumps(res, indent=1, default=str), flush=True)
    res["completion"] = completion()
    res["seconds"] = round(time.time() - t0, 1)
    (OUT / "r07_p22.json").write_text(json.dumps(res, indent=1, default=str) + "\n")
    print(json.dumps(res, indent=1, default=str))


if __name__ == "__main__":
    main()
