#!/usr/bin/env python3
"""XVI-R R04: derive sigma-bar(Y) for Y = c*.G(c,h) from the reduced master S-bar = S_EH + int Ktilde.s0 h + int c*.s0 c,
Ktilde = K + cbar.F'(h).  Three-point pieces:
  X_K : (delta S / delta c)|_K . G = (D^T K).G
  X_C : (delta S / delta c)|_cbar . G = (s F / delta structure).G = ghost kinetic on G
  X_S : - s0 h . (delta Y / delta h)   [h -> s0 h on G, result in c*cc]
The Q08 X_C = (-1/2) p^2 (c*->cbar) is a HAND-MADE kinetic; Q11 replaced it by rho*T(X_K) (ghost-eq imposition).
Here the kinetic coefficient is DERIVED from V_cb / the free identity sF = Box c, and a joint nullspace over
(X_K, X_C_kin, X_S) decides whether a single global (a,b,c) closes every G image on the reduced system."""
import json, pickle, sys, time  # noqa: E401
from pathlib import Path
import sympy as sp
from sympy import QQ_I
from sympy.polys.matrices import DomainMatrix
R = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(R / "r03/code"))
import r03_fix as F  # noqa: E402
S = F.S
I = sp.I
OUT = R / "build" / "logs"
DM, rk = F.DM, F.rk


def pieces():
    """per G monomial: dict of raw vectors X_K, X_C_kin (unit p^2), X_S (h->s0h antisym)."""
    out = []
    for m in S.monomials({"gv": 1, "gc1": 1, "a3": 2}, 2):
        # X_K: (D^T K).G -- same construction as Q11 (Fourier of -partial.K)
        xk = sp.expand(-I * S.U.H("aK", "p2") * S.hhc.subst_vec(m, "gv", {"aK": 1}))
        # X_C kinetic unit: Box on the c* momentum -> (p2^2) (c* -> cbar); coefficient free
        xc1 = sp.expand(S.U.H("p2", "p2") * S.hhc.subst_vec(m, "gv", {"gb2": 1}))
        # X_S: h -> s0 h = D c
        s = S.R.remap(m, {"p1": S.R.P1, "p2": S.R.P3, "gv": "z"})
        s = S.R.strip_t(-I * S.R.polarise(sp.expand(s), "a3", "p2", "gc2"))
        s_sw = S.R.remap(s, {"p1": S.R.P2, "p2": S.R.P1, "gc1": "gc2", "gc2": "gc1"})
        xs = sp.expand(s - s_sw)
        vk, l1 = S.vec("K", xk)
        vc, l2 = S.vec("C", xc1)
        vs, l3 = S.vec("S", xs)
        out.append({"K": vk, "C": vc, "S": vs, "leak": l1 or l2 or l3,
                    "zero": not any(sp.expand(x) != 0 for x in vk + vc + vs)})
    return out


def main():
    t0 = time.time()
    S.SIG = sp.sympify(json.loads((S.LOG / "q07_hhc_tree.json").read_text())["sigma_F"])
    d = pickle.loads((S.LOG / "q08_triv.pkl").read_bytes())
    B = pickle.loads((S.LOG / "q11_blocks.pkl").read_bytes())
    Ab = DM(pickle.loads((OUT / "r03_reduced_hhc.pkl").read_bytes()))
    AAb = Ab.vstack(DM(B["ccK"])).vstack(DM(B["ccb"])).vstack(DM(B["bose"])).vstack(DM(d["Ag"]))
    ps = pieces()
    res = {"n": len(ps), "leak": sum(p["leak"] for p in ps), "n_nonzero": sum(not p["zero"] for p in ps)}
    # joint nullspace: AAb * (a X_K + b X_C + c X_S) = 0 for all generators
    cols = []
    for name in ("K", "C", "S"):
        M = (AAb * DM(sp.Matrix([p[name] for p in ps]).T)).to_Matrix()
        cols.append(M.reshape(M.rows * M.cols, 1))
    ns = DM(sp.Matrix.hstack(*cols)).nullspace().to_Matrix()
    res["coeff_nullspace"] = [[str(x) for x in ns.row(i)] for i in range(ns.rows)]
    # known candidates
    cands = {"Q08_handmade": (1, sp.Rational(-1, 2), 0), "Q11_rhoT_style_no_S": None,
             "kinetic_minus1_no_S": (1, -1, 0), "kinetic_plus1_no_S": (1, 1, 0),
             "kinetic_minus1_S_minus1": (1, -1, -1), "kinetic_minus1_S_plus1": (1, -1, 1),
             "kinetic_half_S_minus1": (1, sp.Rational(-1, 2), -1), "kinetic_half_S_plus1": (1, sp.Rational(-1, 2), 1),
             "only_S_minus1": (0, 0, -1), "K_only": (1, 0, 0)}
    # Q11-style: X_C = rho T(X_K), X_S with eps
    q11 = []
    for p, g in zip(ps, S.family_G()):
        xk = sp.Matrix(p["K"])
        # rebuild rho T partner from the expression level: use family_G anti C part
        a = g["anti"]
        xc = [a[i] if S.OFF["C"] <= i < S.OFF["S"] else 0 for i in range(S.NDIM)]
        q11.append((p["K"], xc, g["field"]))
    for eps in (-1, 1):
        vs = [S.add(k, S.add(c, S.scale(eps, s))) for k, c, s in q11]
        P = AAb * DM(sp.Matrix(vs).T)
        cands[f"Q11_rhoT_S_eps{eps:+d}"] = ("eval", sum(not P.extract(list(range(P.shape[0])), [j]).is_zero_matrix
                                                         for j in range(P.shape[1])))
    def eval_abc(a, b, c):
        vs = [S.add(S.scale(a, p["K"]), S.add(S.scale(b, p["C"]), S.scale(c, p["S"]))) for p in ps]
        vs = [v for v in vs if any(sp.expand(x) != 0 for x in v)]
        P = AAb * DM(sp.Matrix(vs).T)
        return len(vs), sum(not P.extract(list(range(P.shape[0])), [j]).is_zero_matrix for j in range(P.shape[1]))
    tab = {}
    for name, coef in cands.items():
        if isinstance(coef, tuple) and coef[0] == "eval":
            tab[name] = {"failing": coef[1]}
        elif coef is None:
            continue
        else:
            n, f = eval_abc(*coef)
            tab[name] = {"a_b_c": [str(x) for x in coef], "n_images": n, "failing": f}
    res["candidates"] = tab
    # if nullspace non-empty, evaluate each basis vector
    for i in range(ns.rows):
        a, b, c = (sp.sympify(ns[i, j]) for j in range(3))
        n, f = eval_abc(a, b, c)
        res[f"null_{i}"] = {"a_b_c": [str(a), str(b), str(c)], "n_images": n, "failing": f}
    res["seconds"] = round(time.time() - t0, 1)
    (OUT / "r04_gderive.json").write_text(json.dumps(res, indent=1) + "\n")
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
