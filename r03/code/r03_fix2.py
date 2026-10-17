#!/usr/bin/env python3
"""XVI-R R03 step 2 (reduced formulation, continued):
 (a) block-resolve the sigma-bar failures (reduced hhc, ccK, ccb, Bose, ghost eq);
 (b) is there ONE scalar lam with (gf-channel + lam * E-channel) = 0 on ker Ag (i.e. full row == reduced row under the
     ghost equation)?  lam = sigma_F would mean the frozen rows are consistent;
 (c) the gh -1 three-point family Ktilde.Ktilde.c (missing from Q08..R02): sigma-bar replaces one Ktilde by
     dS_EH/dh -> X_Khc, with cbar partner rho T(X_K); f antisymmetric (two Grassmann-odd K's);
 (d) reduced cohomology: Nbar_g / (sigma-bar span cap Nbar_g + blind-bar_g), classification of the 43 transported."""
import json, pickle, sys, time  # noqa: E401
from pathlib import Path
import sympy as sp
from sympy import QQ_I
from sympy.polys.matrices import DomainMatrix
R = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import r03_fix as F  # noqa: E402
S = F.S
I = sp.I
OUT = R / "build" / "logs"
DM, rk, hcat = F.DM, F.rk, F.hcat


def blocks_gf_E():
    XK, xk = F.ansatz({"aK": 2, "gc1": 1, "a3": 2}, 3, "xk")
    XC, xc = F.ansatz({"gb2": 1, "gc1": 1, "a3": 2}, 4, "xc")
    XS, xs = F.ansatz({"z": 1, "gc1": 1, "gc2": 1}, 3, "xs")
    XH, xh = F.ansatz({"a1": 2, "a2": 2, "a3": 2}, 4, "xh")
    xall = xk + xc + xs + xh
    rg, re = 0, 0
    for i in ("2", "3"):
        ai, ki, Fi = S.hhc.EXT[i]
        Gk = S.R.polarise(sp.expand(S.hhc._relabel(XK, i)), "aK", "KU", "KW")
        eh, gf = S.G.S2_EH_gf(S.hhc.ext(ai), S.G.probe_mat("E"), ki)
        rg += S.R.strip_t(S.R.contract_id(Gk, I * gf, "K", "E", "h"))
        re += S.R.strip_t(S.hhc.subst_vec(sp.expand(S.hhc._relabel(XC, i)), "gb2", S.G.Dvec(S.hhc.ext(ai), Fi)))
    rg, re = sp.expand(rg), sp.expand(re)
    X = set(xall)
    gens = sorted({s_ for r in (rg, re) for s_ in r.free_symbols if s_ not in X}, key=str)
    tot = sp.Poly(rg + sp.Symbol("LAMX") * re, *gens, sp.Symbol("LAMX"))
    Mg = sp.linear_eq_to_matrix(sp.Poly(rg, *gens).coeffs(), xall)[0] if rg != 0 else None
    # common generator monomials for both: use coefficient extraction on a shared monomial list
    monos = sorted(set(sp.Poly(rg, *gens).monoms()) | set(sp.Poly(re, *gens).monoms()))
    def mat(r):
        P = sp.Poly(r, *gens)
        d = dict(zip(P.monoms(), P.coeffs()))
        return sp.linear_eq_to_matrix([d.get(m, 0) for m in monos], xall)[0]
    del tot, Mg
    return mat(rg), mat(re)


def family_KKc():
    out = []
    sw = {"p1": S.R.P1, "p2": S.R.P3, "aK": "a3", "a3": "aK"}
    seen = set()
    for m in S.monomials({"aK": 2, "a3": 2, "gc1": 1}, 1):
        f = sp.expand(m - S.R.remap(m, sw))
        if f == 0 or f in seen or -f in seen:
            continue
        seen.add(f)
        fp = S.R.polarise(f, "a3", "KU", "KW")
        eh, gf = S.G.S2_EH_gf(S.hhc.ext("a3"), S.G.probe_mat("E"), S.hhc.p3)
        xk = S.R.strip_t(S.R.contract_id(fp, I * eh, "K", "E", "h"))
        vk, l1 = S.vec("K", xk)
        vc, l2 = S.vec("C", F.partner_C(xk))
        out.append({"v": S.add(vk, vc), "leak": l1 or l2})
    return out


def main():
    t0 = time.time()
    S.SIG = sp.sympify(json.loads((S.LOG / "q07_hhc_tree.json").read_text())["sigma_F"])
    d = pickle.loads((S.LOG / "q08_triv.pkl").read_bytes())
    B = pickle.loads((S.LOG / "q11_blocks.pkl").read_bytes())
    A, Ag = DM(d["A"]), DM(d["Ag"])
    n = A.shape[1]
    Ab = DM(pickle.loads((OUT / "r03_reduced_hhc.pkl").read_bytes()))
    blk = {"hhc_reduced": Ab, "ccK": DM(B["ccK"]), "ccb": DM(B["ccb"]), "bose": DM(B["bose"]), "ghost_eq": Ag}
    sb = pickle.loads((OUT / "r03_sigma_bar.pkl").read_bytes())
    cols, tags = list(sb["cols"]), list(sb["tags"])
    res = {}
    # (a)
    fails = {}
    for j, v in enumerate(cols):
        dv = DM(sp.Matrix(v))
        k = " ".join(b for b, M in blk.items() if not (M * dv).is_zero_matrix)
        if k:
            fails[f"{tags[j]}: {k}"] = fails.get(f"{tags[j]}: {k}", 0) + 1
    res["sigma_bar_fail_blocks"] = fails
    # (b)
    Mg, Me = blocks_gf_E()
    kg = Ag.nullspace().transpose()
    G1, E1 = (DM(Mg) * kg).to_Matrix(), (DM(Me) * kg).to_Matrix()
    lam = None
    for i in range(E1.rows):
        for j in range(E1.cols):
            if E1[i, j] != 0:
                lam = sp.nsimplify(-G1[i, j] / E1[i, j])
                break
        if lam is not None:
            break
    res["lambda_candidate"] = str(lam)
    res["sigma_F"] = str(S.SIG)
    res["lambda_closes_on_kerAg"] = bool(lam is not None and (G1 + lam * E1).applyfunc(sp.expand).is_zero_matrix)
    res["rank_gf_on_kerAg"], res["rank_E_on_kerAg"] = DM(G1).rank(), DM(E1).rank()
    # (c)
    kkc = family_KKc()
    res["KKc"] = {"n": len(kkc), "leak": sum(g["leak"] for g in kkc)}
    kf = {}
    for g in kkc:
        dv = DM(sp.Matrix(g["v"]))
        if dv.is_zero_matrix:
            continue
        k = " ".join(b for b, M in blk.items() if not (M * dv).is_zero_matrix) or "pass"
        kf[k] = kf.get(k, 0) + 1
    res["KKc"]["blocks"] = kf
    cols += [g["v"] for g in kkc if any(x != 0 for x in g["v"])]
    tags += ["KKc"] * (len(cols) - len(tags))
    # (d) reduced cohomology
    Abar = Ab.vstack(blk["ccK"]).vstack(blk["ccb"]).vstack(blk["bose"])
    AAb = Abar.vstack(Ag)
    Nb = AAb.nullspace().transpose()
    Sig = DM(sp.Matrix(cols).T)
    Z = (AAb * Sig).nullspace().transpose()
    T = Sig * Z
    rows_all = Abar
    zc = [j for j in range(n) if rows_all.extract(list(range(rows_all.shape[0])), [j]).is_zero_matrix]
    Bm = DM(sp.Matrix([[1 if i == j else 0 for j in zc] for i in range(n)]))
    kb = (Ag * Bm).nullspace().transpose()
    Bg = Bm * kb
    dNb = Nb.shape[1]
    dW = rk(hcat(T, Bg))
    res["reduced"] = {"dim_Nbar_g": dNb, "n_blind_coords": len(zc), "dim_blind_in_Nbar_g": rk(Bg),
                      "rank_Sigma_bar": rk(Sig), "dim_T_bar": rk(T), "sigma_bar_inside_Nbar": rk(Sig) == rk(T),
                      "dim_trivial_plus_blind": dW, "dim_H_bar_candidates": dNb - dW,
                      "dim_Nbar_cap_span_Sigma_plus_blind": rk(hcat(Sig, Bg)) + dNb - rk(hcat(Sig, Bg, Nb))}
    # transport of Q08 N_g and its 43
    Ng = A.vstack(Ag).nullspace().transpose()
    res["reduced"]["dim_N_g_cap_Nbar_g"] = Ng.shape[1] + dNb - rk(hcat(Ng, Nb))
    pickle.dump({"cols": cols, "tags": tags}, open(OUT / "r03_sigma_bar2.pkl", "wb"))
    res["seconds"] = round(time.time() - t0, 1)
    (OUT / "r03_fix2.json").write_text(json.dumps(res, indent=1, default=str) + "\n")
    print(json.dumps(res, indent=1, default=str))


if __name__ == "__main__":
    main()
