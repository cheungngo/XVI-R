#!/usr/bin/env python3
"""XVI-R R03: reduced-functional (Gamma-bar) formulation.  With B integrated out the full sigma = (S0, .) is nilpotent only
up to antighost equations of motion; the leg analysis puts the whole hhc defect in A-channel(D_gf^(1) + E^(1)).
Standard BV: Gamma = Gamma-bar + S_gf, and under the ghost equation cbar enters Gamma-bar only through Ktilde, so the
ST identity is ST-bar(Gamma-bar), whose hhc row has channels A and D_EH only (no gauge-fixing kinetic term, no E channel),
and the trivial counterterms are sigma-bar(Y) = (S_EH + antifield terms, Y) on Y = Y(Ktilde, c*, h, c).
Steps (no fitting):
 (1) build the reduced hhc row block A-bar_hhc;
 (2) EQUIVALENCE: on ker(Ag) the full hhc block equals A-bar_hhc (rank of difference restricted to ker Ag);
 (3) sigma-bar images: F: X_H = sum_legs D_EH, X_K = A_K, X_C = rho T(X_K);  G: unchanged (no D/E pieces);
     P (standalone cbar) has no reduced counterpart and is dropped;
 (4) test images in N_g (full A + Ag) and in N-bar_g (reduced), then recompute trivial span and the 43 classification."""
import json, pickle, sys, time  # noqa: E401
from pathlib import Path
import sympy as sp
from sympy import QQ_I
from sympy.polys.matrices import DomainMatrix
R = Path(__file__).resolve().parents[2]
QQ = R / "build" / "xviq" / "XVI-Q"
sys.path[:0] = [str(QQ / "code/p6q06/q07/q08/q09/q10/q11"), str(QQ / "code/p6q06/q07")]
import sigma11 as S  # noqa: E402
from evfit import ansatz  # noqa: E402
I = sp.I
OUT = R / "build" / "logs"


def DM(M):
    return DomainMatrix.from_Matrix(M).convert_to(QQ_I)


def rk(D):
    return D.rank() if D.shape[0] and D.shape[1] else 0


def hcat(*Ds):
    Ds = [d for d in Ds if d.shape[1]]
    return Ds[0].hstack(*Ds[1:]) if len(Ds) > 1 else Ds[0]


def reduced_hhc_block():
    XK, xk = ansatz({"aK": 2, "gc1": 1, "a3": 2}, 3, "xk")
    XC, xc = ansatz({"gb2": 1, "gc1": 1, "a3": 2}, 4, "xc")
    XS, xs = ansatz({"z": 1, "gc1": 1, "gc2": 1}, 3, "xs")
    XH, xh = ansatz({"a1": 2, "a2": 2, "a3": 2}, 4, "xh")
    xall = xk + xc + xs + xh
    r = S.R.strip_t(-I * S.R.polarise(sp.expand(XH), "a1", "p1", "gc1"))
    for i in ("2", "3"):
        ai, ki, Fi = S.hhc.EXT[i]
        Gk = S.R.polarise(sp.expand(S.hhc._relabel(XK, i)), "aK", "KU", "KW")
        eh, gf = S.G.S2_EH_gf(S.hhc.ext(ai), S.G.probe_mat("E"), ki)
        r += S.R.strip_t(S.R.contract_id(Gk, I * eh, "K", "E", "h"))
    r = sp.expand(r)
    X = set(xall)
    gens = sorted([s for s in r.free_symbols if s not in X], key=str)
    return sp.linear_eq_to_matrix(sp.Poly(r, *gens).coeffs(), xall)[0]


def diff_block():
    """full minus reduced hhc row = gauge-fixing D channel on X_K + sigma_F * E channel on X_C (legs 2,3)."""
    XK, xk = ansatz({"aK": 2, "gc1": 1, "a3": 2}, 3, "xk")
    XC, xc = ansatz({"gb2": 1, "gc1": 1, "a3": 2}, 4, "xc")
    XS, xs = ansatz({"z": 1, "gc1": 1, "gc2": 1}, 3, "xs")
    XH, xh = ansatz({"a1": 2, "a2": 2, "a3": 2}, 4, "xh")
    xall = xk + xc + xs + xh
    r = 0
    for i in ("2", "3"):
        ai, ki, Fi = S.hhc.EXT[i]
        Gk = S.R.polarise(sp.expand(S.hhc._relabel(XK, i)), "aK", "KU", "KW")
        eh, gf = S.G.S2_EH_gf(S.hhc.ext(ai), S.G.probe_mat("E"), ki)
        r += S.R.strip_t(S.R.contract_id(Gk, I * gf, "K", "E", "h"))
        r += S.SIG * S.R.strip_t(S.hhc.subst_vec(sp.expand(S.hhc._relabel(XC, i)), "gb2", S.G.Dvec(S.hhc.ext(ai), Fi)))
    r = sp.expand(r)
    X = set(xall)
    gens = sorted([s_ for s_ in r.free_symbols if s_ not in X], key=str)
    return sp.linear_eq_to_matrix(sp.Poly(r, *gens).coeffs(), xall)[0]


def partner_C(xk):
    xk = sp.nsimplify(sp.expand(xk), rational=True)
    return sp.expand(S.RHO * S.Tmap(xk)) if xk != 0 else sp.Integer(0)


def main():
    t0 = time.time()
    S.SIG = sp.sympify(json.loads((S.LOG / "q07_hhc_tree.json").read_text())["sigma_F"])
    d = pickle.loads((S.LOG / "q08_triv.pkl").read_bytes())
    B = pickle.loads((S.LOG / "q11_blocks.pkl").read_bytes())
    A, Ag = DM(d["A"]), DM(d["Ag"])
    n = A.shape[1]
    Ab = reduced_hhc_block()
    pickle.dump(Ab, open(OUT / "r03_reduced_hhc.pkl", "wb"))
    Dif = DM(diff_block())
    Ab = DM(Ab)
    others = [DM(B[k]) for k in ("ccK", "ccb", "bose")]
    Abar_full = Ab
    for o in others:
        Abar_full = Abar_full.vstack(o)
    res = {}
    kg = Ag.nullspace().transpose()
    res["equivalence_rank_full_minus_reduced_on_kerAg"] = rk(Dif * kg)
    res["dim_kerAg"] = kg.shape[1]
    print(res, flush=True)
    Ng = A.vstack(Ag).nullspace().transpose()
    Nbg = Abar_full.vstack(Ag).nullspace().transpose()
    res["dim_N_g"], res["dim_Nbar_g"] = Ng.shape[1], Nbg.shape[1]
    res["dim_N_g_plus_Nbar_g"] = rk(hcat(Ng, Nbg))
    # sigma-bar images
    cols, tags = [], []
    for m in S.monomials({"aK": 2, "a2": 2, "a3": 2}, 2):
        f = sp.expand(m + S.R.remap(m, S.SW23))
        if f == 0:
            continue
        xh = 0
        for i in "123":
            fi = S.R.remap(f, S.REL[i]) if S.REL[i] else f
            fi = S.R.polarise(fi, "aK", "KU", "KW")
            eh, gf = S.G.S2_EH_gf(S.hhc.ext(f"a{i}"), S.G.probe_mat("E"), S.MOM[i])
            xh += S.R.strip_t(S.R.contract_id(fi, I * eh, "K", "E", "h"))
        xk = S.A_on_K(f)
        v = S.add(S.vec("H", xh)[0], S.scale(-1, S.add(S.vec("K", xk)[0], S.vec("C", partner_C(xk))[0])))
        if any(sp.expand(x) != 0 for x in v):
            cols.append(v)
            tags.append("F")
    for g in S.family_G():
        v = S.add(g["anti"], S.scale(-1, g["field"]))
        if any(sp.expand(x) != 0 for x in v):
            cols.append(v)
            tags.append("G")
    Sig = DM(sp.Matrix(cols).T)
    AAf, AAb = A.vstack(Ag), Abar_full.vstack(Ag)
    for lab, M in (("full", AAf), ("reduced", AAb)):
        P = M * Sig
        bad = [j for j in range(P.shape[1]) if not P.extract(list(range(P.shape[0])), [j]).is_zero_matrix]
        res[f"images_failing_{lab}"] = {"F": sum(tags[j] == "F" for j in bad), "G": sum(tags[j] == "G" for j in bad)}
        res[f"defect_rank_{lab}"] = rk(P)
    res["n_images"] = {"F": tags.count("F"), "G": tags.count("G")}
    res["rank_Sigma_bar"] = rk(Sig)
    pickle.dump({"cols": cols, "tags": tags}, open(OUT / "r03_sigma_bar.pkl", "wb"))
    # classification of the 43 with sigma-bar (constraints: full A + Ag, i.e. N_g of Q08)
    zc = [j for j in range(n) if all(d["A"][i, j] == 0 for i in range(d["A"].shape[0]))]
    Bm = DM(sp.Matrix([[1 if i == j else 0 for j in zc] for i in range(n)]))
    Bg = Bm * (Ag * Bm).nullspace().transpose()
    keep = pickle.loads((S.LOG / "q11_sigma_vecs.pkl").read_bytes())
    W11 = hcat(DM(sp.Matrix([*keep[("F", "derived_minus1")], *keep[("G", "derived_minus1")]]).T), Bg)
    Z = (AAf * Sig).nullspace().transpose()
    T = Sig * Z if Z.shape[1] else DomainMatrix.zeros((n, 0), QQ_I)
    res["dim_T_bar_in_Ng"] = rk(T)
    res["Sigma_bar_inside_Ng"] = rk(Sig) == rk(T)
    # R02 trivial span (computed sigma incl P) for comparison
    r02 = pickle.loads((OUT / "r02_T.pkl").read_bytes()) if (OUT / "r02_T.pkl").is_file() else None
    Wp = hcat(W11, T, Bg) if T.shape[1] else hcat(W11, Bg)
    dNg, dW11, dWp = Ng.shape[1], rk(W11), rk(hcat(W11, T, Bg) if T.shape[1] else hcat(W11, Bg))
    res["dim_W_Q11"], res["dim_Wbar_cap_Ng"] = dW11, dWp
    res["dim_Ng_cap_span_Sigma_bar_plus_blind"] = rk(hcat(Sig, Bg)) + dNg - rk(hcat(Sig, Bg, Ng))
    res["classification_43_sigma_bar"] = {"trivial": dWp - dW11, "blind": 0, "candidates": dNg - dWp}
    res["seconds"] = round(time.time() - t0, 1)
    (OUT / "r03_fix.json").write_text(json.dumps(res, indent=1) + "\n")
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
