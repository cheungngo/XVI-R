#!/usr/bin/env python3
"""XVI-R R04 final: derived sigma-bar on Y=c*.G from the reduced master; G-sector outcome; classification.
Derived (nothing fitted):
  S-bar = S_EH + int Ktilde.s0 h + int c*.s0 c,  Ktilde = K + cbar F'(h)  (cbar not independent).
  For Y = c*.G(c,h):
    X_K  = (D^T Ktilde).G     [delta S / delta c |_Ktilde . G]
    X_C  = rho * T(X_K)         [expansion of Ktilde in (K,cbar); NOT an independent kinetic]
    X_S  = - s0 h . dY/dh     [h -> s0 h; eps_S = -1 from the antibracket field/antifield relative sign]
  The raw kinetic alternative cbar.Box.G (Q08 handmade) fails all 33; coefficient scan on (X_K, kinetic, X_S) is empty.
  With the derived image: 20/33 pass the reduced system+ghost eq; 13/33 fail both ghost rows ccK and ccb.
  Pattern: every failure has a single hatted contraction on the c* index (T_*gv with only H elsewhere);
  the 6 double-T monomials with T_*gv pass.  No completion of the three pieces removes the 13.
  Obstruction: OPEN (G-sector T(c*)), conditional classes remain conditional."""
import json, pickle, sys, time  # noqa: E401
from pathlib import Path
import sympy as sp
from sympy import QQ_I
from sympy.polys.matrices import DomainMatrix
R = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(R / "r03/code"))
import r03_fix as F  # noqa: E402
sys.path.insert(0, str(R / "code/r02"))
from r02_sci import tex_mono  # noqa: E402
S = F.S
OUT = R / "build" / "logs"
DM, rk, hcat = F.DM, F.rk, F.hcat


def coef(c):
    import re
    c = c.replace(" ", "")
    c = re.sub(r"\\frac\{([^}]*)\}\{([^}]*)\}", r"(\1)/(\2)", c)
    c = re.sub(r"(\d)i", r"\1*I", c).replace("(i)", "(I)")
    c = "I" if c == "i" else ("-I" if c == "-i" else c)
    return sp.sympify(c, locals={"I": sp.I})


def G_image(m, eps_S=-1):
    xk = sp.expand(-sp.I * S.U.H("aK", "p2") * S.hhc.subst_vec(m, "gv", {"aK": 1}))
    xc = sp.expand(S.RHO * S.Tmap(sp.nsimplify(xk, rational=True))) if xk != 0 else 0
    s = S.R.remap(m, {"p1": S.R.P1, "p2": S.R.P3, "gv": "z"})
    s = S.R.strip_t(-sp.I * S.R.polarise(sp.expand(s), "a3", "p2", "gc2"))
    xs = sp.expand(s - S.R.remap(s, {"p1": S.R.P2, "p2": S.R.P1, "gc1": "gc2", "gc2": "gc1"}))
    return S.add(S.vec("K", xk)[0], S.vec("C", xc)[0] if xc != 0 else [0] * S.NDIM,
                 S.scale(eps_S, S.vec("S", xs)[0]))


def main():
    t0 = time.time()
    S.SIG = sp.sympify(json.loads((S.LOG / "q07_hhc_tree.json").read_text())["sigma_F"])
    d = pickle.loads((S.LOG / "q08_triv.pkl").read_bytes())
    B = pickle.loads((S.LOG / "q11_blocks.pkl").read_bytes())
    A, Ag = DM(d["A"]), DM(d["Ag"])
    n = A.shape[1]
    Ab = DM(pickle.loads((OUT / "r03_reduced_hhc.pkl").read_bytes()))
    Abar = Ab.vstack(DM(B["ccK"])).vstack(DM(B["ccb"])).vstack(DM(B["bose"]))
    AAb = Abar.vstack(Ag)
    sb = pickle.loads((OUT / "r03_sigma_bar2.pkl").read_bytes())
    Fc = [v for v, t in zip(sb["cols"], sb["tags"]) if t == "F"]
    Kc = [v for v, t in zip(sb["cols"], sb["tags"]) if t == "KKc"]
    # G derived
    gd = json.loads((OUT / "r04_gderive.json").read_text())
    gp = json.loads((OUT / "r04_gpat.json").read_text())
    Gs, tags_g, mono = [], [], []
    for m in S.monomials({"gv": 1, "gc1": 1, "a3": 2}, 2):
        v = G_image(m, -1)
        if not any(sp.expand(x) != 0 for x in v):
            continue
        ok = (AAb * DM(sp.Matrix(v))).is_zero_matrix
        Gs.append(v)
        tags_g.append("pass" if ok else "fail")
        mono.append(str(m))
    res = {"derivation": {"X_K": "(D^T Ktilde).G", "X_C": "rho*T(X_K) (Ktilde expansion; not independent kinetic)",
                          "X_S": "-s0h.dY/dh with eps_S=-1 (antibracket relative sign)",
                          "kinetic_alternative": "cbar.Box.G fails all 33; nullspace(X_K,kinetic,X_S) empty",
                          "eps_S": -1},
           "gderive_nullspace": gd["coeff_nullspace"], "gderive_candidates": gd["candidates"],
           "G": {"n_images": len(Gs), "n_pass": tags_g.count("pass"), "n_fail": tags_g.count("fail"),
                 "fail_blocks": "ccK+ccb for all failures", "pattern": gp["by_fail"],
                 "failures": gp["failures"],
                 "outcome": "PARTIAL: 20/33 closed; 13/33 OPEN (T(c*) obstruction)"}}
    # full sigma-bar with only passing G (honest: failing Y excluded from the proven trivial span)
    cols = Fc + [v for v, t in zip(Gs, tags_g) if t == "pass"] + Kc
    tags = (["F"] * len(Fc) + ["G_pass"] * tags_g.count("pass") + ["KKc"] * len(Kc))
    # also span of ALL G (incl fail) for robustness
    cols_all = Fc + Gs + Kc
    Sig = DM(sp.Matrix(cols).T)
    SigA = DM(sp.Matrix(cols_all).T)
    Nb = AAb.nullspace().transpose()
    Z = (AAb * Sig).nullspace().transpose()
    T = Sig * Z
    zc = [j for j in range(n) if Abar.extract(list(range(Abar.shape[0])), [j]).is_zero_matrix]
    Bm = DM(sp.Matrix([[1 if i == j else 0 for j in zc] for i in range(n)]))
    Bg = Bm * (Ag * Bm).nullspace().transpose()
    dNb, dT, dB, dTB = Nb.shape[1], rk(T), rk(Bg), rk(hcat(T, Bg))
    dRob = rk(hcat(SigA, Bg)) + dNb - rk(hcat(SigA, Bg, Nb))
    res["reduced"] = {"dim_Nbar_g": dNb, "rank_Sigma_bar_passing": rk(Sig), "dim_T_bar": dT, "dim_blind_g": dB,
                      "dim_T_cap_blind": dT + dB - dTB, "dim_H_bar": dNb - dTB,
                      "dim_span_allG_plus_blind_cap_Nbar": dRob,
                      "sigma_bar_nilpotent_on_passing_Y": True,
                      "n_Y_passing": len(cols), "n_Y_all": len(cols_all)}
    # canonical basis
    W = hcat(T, Bg)
    Wr, Pw = W.transpose().to_field().rref()
    Wr = Wr.to_Matrix()[:len(Pw), :]
    Nr = Nb.transpose().to_field().to_Matrix()
    for i, pc in enumerate(Pw):
        for r in range(Nr.shape[0]):
            if Nr[r, pc] != 0:
                Nr[r, :] = Nr[r, :] - Nr[r, pc] * Wr[i, :]
    Hr, Ph = DM(Nr).to_field().rref()
    Hr = Hr.to_Matrix()[:len(Ph), :]
    lab = [(sec, m) for sec in "KCSH" for m in S.MS[sec]]
    basis = []
    for r in range(Hr.shape[0]):
        v = Hr[r, :]
        nz = [i for i in range(n) if v[i] != 0]
        terms = [(lab[i][0], tex_mono(lab[i][1]), sp.latex(sp.nsimplify(v[i]))) for i in nz]
        basis.append({"pivot": [lab[Ph[r]][0], str(lab[Ph[r]][1])], "sectors": sorted({t[0] for t in terms}),
                      "n_terms": len(terms), "terms": terms})
    res["H_bar_basis"] = basis
    # transport R03's 18 (= R02's 22 after R03)
    r03b = json.loads((R / "certificates/R03_H_bar_basis.json").read_text())["basis"]
    idx = {(sec, tex_mono(m)): i for i, (sec, m) in enumerate(lab)}
    reps = []
    for b in r03b:
        v = [0] * n
        for sec, m, c in b["terms"]:
            v[idx[(sec, m)]] = coef(c)
        reps.append(v)
    C = DM(sp.Matrix(reps).T)
    dC = rk(C)
    # intersect with Nbar, then with T+blind
    M = hcat(C, Nb)
    ker = M.nullspace().transpose().to_Matrix()
    CNv = (C.to_Matrix() * ker[:C.shape[1], :]) if ker.shape[1] else sp.zeros(n, 0)
    CNd = DM(CNv) if CNv.shape[1] else DomainMatrix.zeros((n, 0), QQ_I)
    triv = (rk(CNd) + dTB - rk(hcat(CNd, W))) if CNd.shape[1] else 0
    res["transport_18"] = {"dim": dC, "inside_Nbar_g": rk(CNd),
                           "trivial_or_blind": triv, "nontrivial_in_H_bar": rk(CNd) - triv,
                           "outside_all_computed_images": rk(CNd) - triv}  # equals H_bar dim if they span it
    res["seconds"] = round(time.time() - t0, 1)
    (OUT / "r04_final.json").write_text(json.dumps(res, indent=1, default=str) + "\n")
    print(json.dumps({k: v for k, v in res.items() if k != "H_bar_basis"}, indent=1, default=str))


if __name__ == "__main__":
    main()
