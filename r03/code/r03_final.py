#!/usr/bin/env python3
"""XVI-R R03 final: reduced (Gamma-bar) cohomology at kappa^3, three points, pole level.
G family: the relative sign eps_S of the c*cc piece was left undetermined in Q11 ('both signs reported'); both are tested
here against the reduced system (sigma-bar^2 = 0 must hold for every generator; one binary choice vs 33 generators).
Y-space: F (Ktilde h h), G (c* c h), KKc (Ktilde Ktilde c) -- every gh -1 three-field content in (Ktilde, c*, h, c).
Outputs: reduced classification and the transport of R02's 22 candidates (canonical basis, in N_g)."""
import json, pickle, sys, time  # noqa: E401
from pathlib import Path
import sympy as sp
from sympy import QQ_I
from sympy.polys.matrices import DomainMatrix
R = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import r03_fix as F  # noqa: E402
import r03_fix2 as F2  # noqa: E402,F401
sys.path.insert(0, str(R / 'code/r02'))
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
    res = {"n_F": len(Fc), "n_KKc": len(Kc)}
    fam = S.family_G()
    Gc = {}
    for eps in (-1, 1):
        vs = [S.add(g["anti"], S.scale(eps, g["field"])) for g in fam]
        vs = [v for v in vs if any(x != 0 for x in v)]
        P = AAb * DM(sp.Matrix(vs).T)
        res[f"G_eps_{eps:+d}_failing"] = sum(not P.extract(list(range(P.shape[0])), [j]).is_zero_matrix for j in range(P.shape[1]))
        Gc[eps] = vs
    eps = min((-1, 1), key=lambda e: res[f"G_eps_{e:+d}_failing"])
    res["eps_S_selected"] = eps
    res["eps_S_unique"] = res["G_eps_-1_failing"] != res["G_eps_+1_failing"]
    cols = Fc + Gc[eps] + Kc
    Sig = DM(sp.Matrix(cols).T)
    P = AAb * Sig
    res["n_images_failing_reduced"] = sum(not P.extract(list(range(P.shape[0])), [j]).is_zero_matrix for j in range(P.shape[1]))
    res["sigma_bar_nilpotent_on_Y"] = res["n_images_failing_reduced"] == 0
    Nb = AAb.nullspace().transpose()
    Z = (AAb * Sig).nullspace().transpose()
    T = Sig * Z
    zc = [j for j in range(n) if Abar.extract(list(range(Abar.shape[0])), [j]).is_zero_matrix]
    Bm = DM(sp.Matrix([[1 if i == j else 0 for j in zc] for i in range(n)]))
    Bg = Bm * (Ag * Bm).nullspace().transpose()
    dNb, dT, dB, dTB = Nb.shape[1], rk(T), rk(Bg), rk(hcat(T, Bg))
    res["reduced"] = {"dim_Nbar_g": dNb, "rank_Sigma_bar": rk(Sig), "dim_T_bar": dT, "dim_blind_g": dB,
                      "dim_T_cap_blind": dT + dB - dTB, "dim_H_bar_nonblind": dNb - dTB,
                      "dim_span_Sigma_plus_blind_cap_Nbar": rk(hcat(Sig, Bg)) + dNb - rk(hcat(Sig, Bg, Nb))}
    # canonical basis of H-bar (rref of Nbar modulo T+blind)
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
    lab = []
    for sec in "KCSH":
        lab += [(sec, m) for m in S.MS[sec]]
    basis = []
    for r in range(Hr.shape[0]):
        v = Hr[r, :]
        nz = [i for i in range(n) if v[i] != 0]
        terms = [(lab[i][0], tex_mono(lab[i][1]), sp.latex(sp.nsimplify(v[i]))) for i in nz]
        basis.append({"pivot": [lab[Ph[r]][0], str(lab[Ph[r]][1])], "sectors": sorted({t[0] for t in terms}),
                      "n_terms": len(terms), "terms": terms})
    res["H_bar_basis"] = basis
    # transport of R02's 22 candidates (vectors from the R02 canonical basis, rebuilt from its certificate order)
    r02 = json.loads((R / "certificates/R02_candidate_basis.json").read_text())
    idx = {(sec, tex_mono(m)): i for i, (sec, m) in enumerate(lab)}
    reps = []
    for b in r02["basis"]:
        v = [0] * n
        for sec, m, c in b["terms"]:
            v[idx[(sec, m)]] = coef(c)
        reps.append(v)
    res["R02_basis_reparsed"] = len(reps)
    if reps is not None:
        C = DM(sp.Matrix(reps).T)
        dC = rk(C)
        dCN = dC + dNb - rk(hcat(C, Nb))
        # dimension of (C cap Nbar) cap (T + blind)
        M = hcat(C, Nb)
        ker = M.nullspace().transpose().to_Matrix()
        CNv = (C.to_Matrix() * ker[:C.shape[1], :]) if ker.shape[1] else sp.zeros(n, 0)
        CNd = DM(CNv) if CNv.shape[1] else None
        triv = (rk(CNd) + dTB - rk(hcat(CNd, W))) if CNd is not None else 0
        res["transport_22"] = {"dim": dC, "inside_Nbar_g": dCN, "fail_reduced_constraints": dC - dCN,
                               "trivial_or_blind_in_reduced": triv, "nontrivial_classes_in_H_bar": dCN - triv}
    res["seconds"] = round(time.time() - t0, 1)
    (OUT / "r03_final.json").write_text(json.dumps(res, indent=1, default=str) + "\n")
    print(json.dumps({k: v for k, v in res.items() if k != "H_bar_basis"}, indent=1, default=str))


if __name__ == "__main__":
    main()
