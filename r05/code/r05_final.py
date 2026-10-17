#!/usr/bin/env python3
"""XVI-R R05 final: reduced system with hat-aware ghost rows (ccK', ccb'); G/F/KKc images; classification; transport
of R04's 20 classes.  Sanity: the non-hat rebuilt blocks reproduce the frozen Q11 blocks (row-space equality)."""
import json, pickle, sys, time  # noqa: E401
from pathlib import Path
import sympy as sp
from sympy import QQ_I
from sympy.polys.matrices import DomainMatrix
R = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(R / "r04/code"))
import r04_final as F4  # noqa: E402
F = F4.F
S = F.S
OUT = R / "build/logs"
DM, rk, hcat = F.DM, F.rk, F.hcat


def same_rowspace(A, B):
    return rk(A) == rk(B) == rk(A.vstack(B))


def classify(AAb, Abar, Ag, cols, n):
    Sig = DM(sp.Matrix(cols).T)
    Nb = AAb.nullspace().transpose()
    Z = (AAb * Sig).nullspace().transpose()
    T = Sig * Z
    zc = [j for j in range(n) if Abar.extract(list(range(Abar.shape[0])), [j]).is_zero_matrix]
    Bm = DM(sp.Matrix([[1 if i == j else 0 for j in zc] for i in range(n)]))
    Bg = Bm * (Ag * Bm).nullspace().transpose()
    dNb, dT, dB, dTB = Nb.shape[1], rk(T), rk(Bg), rk(hcat(T, Bg))
    return {"dim_Nbar_g": dNb, "rank_Sigma": rk(Sig), "dim_T": dT, "dim_blind_g": dB, "dim_T_cap_blind": dT + dB - dTB,
            "dim_H": dNb - dTB, "sigma_inside_Nbar": rk(Sig) == dT}, Nb, T, Bg


def main():
    t0 = time.time()
    S.SIG = sp.sympify(json.loads((S.LOG / "q07_hhc_tree.json").read_text())["sigma_F"])
    d = pickle.loads((S.LOG / "q08_triv.pkl").read_bytes())
    B = pickle.loads((S.LOG / "q11_blocks.pkl").read_bytes())
    NB = pickle.loads((OUT / "r05_blocks.pkl").read_bytes())
    Ag = DM(d["Ag"])
    n = Ag.shape[1]
    res = {"sanity_old_blocks_equal_Q11": {t: same_rowspace(DM(NB[t + "_old"]), DM(B[t])) for t in ("ccK", "ccb")},
           "hat_changes_rowspace": {t: not same_rowspace(DM(NB[t]), DM(B[t])) for t in ("ccK", "ccb")}}
    Ab = DM(pickle.loads((OUT / "r03_reduced_hhc.pkl").read_bytes()))
    Abar = Ab.vstack(DM(NB["ccK"])).vstack(DM(NB["ccb"])).vstack(DM(B["bose"]))
    AAb = Abar.vstack(Ag)
    blk = {"hhc_red": Ab, "ccK_hat": DM(NB["ccK"]), "ccb_hat": DM(NB["ccb"]), "bose": DM(B["bose"]), "geq": Ag}
    sb = pickle.loads((OUT / "r03_sigma_bar2.pkl").read_bytes())
    Fc = [v for v, t in zip(sb["cols"], sb["tags"]) if t == "F"]
    Kc = [v for v, t in zip(sb["cols"], sb["tags"]) if t == "KKc"]
    fam = {"F": Fc, "KKc": Kc, "G": []}
    gm = []
    for m in S.monomials({"gv": 1, "gc1": 1, "a3": 2}, 2):
        v = F4.G_image(m, -1)
        if any(sp.expand(x) != 0 for x in v):
            fam["G"].append(v)
            gm.append(str(m))
    fails = {}
    gfail = []
    for k, vs in fam.items():
        nf = 0
        for i, v in enumerate(vs):
            dv = DM(sp.Matrix(v))
            bad = [b for b, M in blk.items() if not (M * dv).is_zero_matrix]
            if bad:
                nf += 1
                fails[f"{k}: {' '.join(bad)}"] = fails.get(f"{k}: {' '.join(bad)}", 0) + 1
                if k == "G":
                    gfail.append(gm[i])
        res[f"{k}_images"] = {"n": len(vs), "failing": nf}
    res["fail_blocks"] = fails
    res["G_failing_monomials"] = gfail
    # eps_S re-check under the hat-aware rows
    for eps in (1,):
        vs = [F4.G_image(sp.sympify(m), eps) for m in gm]
        P = AAb * DM(sp.Matrix(vs).T)
        res["G_eps_plus1_failing"] = sum(not P.extract(list(range(P.shape[0])), [j]).is_zero_matrix for j in range(P.shape[1]))
    cols = Fc + fam["G"] + Kc
    cl, Nb, T, Bg = classify(AAb, Abar, Ag, cols, n)
    res["reduced_hat"] = cl
    W = hcat(T, Bg)
    res["dim_Nbar_cap_span_all_plus_blind"] = rk(hcat(DM(sp.Matrix(cols).T), Bg)) + Nb.shape[1] - \
        rk(hcat(DM(sp.Matrix(cols).T), Bg, Nb))
    # canonical basis
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
        terms = [(lab[i][0], F4.tex_mono(lab[i][1]), sp.latex(sp.nsimplify(v[i]))) for i in nz]
        basis.append({"pivot": [lab[Ph[r]][0], str(lab[Ph[r]][1])], "sectors": sorted({t[0] for t in terms}),
                      "n_terms": len(terms), "terms": terms})
    res["H_basis"] = basis
    # transport R04's 20
    r04b = json.loads((R / "certificates/R04_H_bar_basis.json").read_text())["basis"]
    idx = {(sec, F4.tex_mono(m)): i for i, (sec, m) in enumerate(lab)}
    reps = []
    for b in r04b:
        v = [0] * n
        for sec, m, c in b["terms"]:
            v[idx[(sec, m)]] = F4.coef(c)
        reps.append(v)
    C = DM(sp.Matrix(reps).T)
    M = hcat(C, Nb)
    ker = M.nullspace().transpose().to_Matrix()
    CNv = (C.to_Matrix() * ker[:C.shape[1], :]) if ker.shape[1] else sp.zeros(n, 0)
    CNd = DM(CNv) if CNv.shape[1] else DomainMatrix.zeros((n, 0), QQ_I)
    dCN = rk(CNd) if CNd.shape[1] else 0
    triv = (dCN + rk(W) - rk(hcat(CNd, W))) if dCN else 0
    res["transport_R04_20"] = {"dim": rk(C), "inside_Nbar_hat": dCN, "not_cocycles": rk(C) - dCN,
                               "trivial_or_blind": triv, "nontrivial": dCN - triv}
    res["seconds"] = round(time.time() - t0, 1)
    (OUT / "r05_final.json").write_text(json.dumps(res, indent=1, default=str) + "\n")
    print(json.dumps({k: v for k, v in res.items() if k != "H_basis"}, indent=1, default=str))


if __name__ == "__main__":
    main()
