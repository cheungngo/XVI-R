#!/usr/bin/env python3
"""XVI-R R03: leg-resolved sub-identities of sigma^2 = 0 on the F family (hhc row).  sigma(K.F) has D-pieces with K at
leg i = 1,2,3 of X_hhh.  In the hhc row the c sits on leg 1 (A channel acts on leg 1 of X_hhh), the D/E channels on legs 2,3.
Expected (nilpotency): (I1) A-channel on D_EH^(1) vanishes (gauge invariance of the EH kinetic term);
(I2) A-channel on D^(2)+D^(3) equals the D channel on A_K (variations on different legs commute);
(I3) the remainder A-channel on D_gf^(1) + E/A_C pieces cancels (gauge-fixing / antighost sector).
Each identity is tested as rank of A applied to the combination over all F generators (0 = holds)."""
import json, pickle, sys, time  # noqa: E401
from pathlib import Path
import sympy as sp
from sympy import QQ_I
from sympy.polys.matrices import DomainMatrix
R = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(R / "build/xviq/XVI-Q/code/p6q06/q07/q08/q09/q10/q11"))
import sigma11 as S  # noqa: E402
I = sp.I


def DM(M):
    return DomainMatrix.from_Matrix(M).convert_to(QQ_I)


def D_leg(f, i, part):
    fi = S.R.remap(f, S.REL[i]) if S.REL[i] else f
    fi = S.R.polarise(fi, "aK", "KU", "KW")
    eh, gf = S.G.S2_EH_gf(S.hhc.ext(f"a{i}"), S.G.probe_mat("E"), S.MOM[i])
    return S.R.strip_t(S.R.contract_id(fi, I * (eh if part == "EH" else gf), "K", "E", "h"))


def E_leg(phi, i):
    fi = S.R.remap(phi, S.REL[i]) if S.REL[i] else phi
    D = S.G.Dvec(S.hhc.ext(f"a{i}"), S.U.V({"p1": 1}) if i == "1" else (S.hhc.F2 if i == "2" else S.hhc.F3))
    return S.R.strip_t(S.SIG * S.hhc.subst_vec(fi, "gv", D))


def main():
    t0 = time.time()
    S.SIG = sp.sympify(json.loads((S.LOG / "q07_hhc_tree.json").read_text())["sigma_F"])
    A = DM(pickle.loads((S.LOG / "q11_blocks.pkl").read_bytes())["hhc"])   # hhc row block only (Bose excluded)
    z = [0] * S.NDIM
    P = {}
    for m in S.monomials({"aK": 2, "a2": 2, "a3": 2}, 2):
        f = sp.expand(m + S.R.remap(m, S.SW23))
        if f == 0:
            continue
        phi = S.cbar_partner(f)
        for i in "123":
            for part in ("EH", "gf"):
                P.setdefault(f"D_{part}{i}", []).append(S.vec("H", D_leg(f, i, part))[0])
            P.setdefault(f"E{i}", []).append(S.vec("H", E_leg(phi, i))[0] if phi != 0 else z)
        P.setdefault("A_K", []).append(S.vec("K", S.A_on_K(f))[0])
        P.setdefault("A_C", []).append(S.vec("C", S.A_on_C(phi))[0] if phi != 0 else z)
    B = {k: A * DM(sp.Matrix(v).T) for k, v in P.items()}

    def rk(*terms):
        M = None
        for c, k in terms:
            t = B[k] * DomainMatrix([[QQ_I.convert(c)]], (1, 1), QQ_I) if False else B[k]
            t = t if c == 1 else (-t if c == -1 else t)
            M = t if M is None else M + t
        return M.rank()
    res = {"rank_each": {k: v.rank() for k, v in B.items()},
           "I1_EH_leg1": rk((1, "D_EH1")),
           "I2_legs23_vs_A_K": rk((1, "D_EH2"), (1, "D_EH3"), (1, "D_gf2"), (1, "D_gf3"), (-1, "A_K")),
           "I2_EH_only_legs23_vs_A_K": rk((1, "D_EH2"), (1, "D_EH3"), (-1, "A_K")),
           "I3_gf1_E_AC": rk((1, "D_gf1"), (1, "E1"), (1, "E2"), (1, "E3"), (-1, "A_C")),
           "I3_E_legs23_vs_A_C": rk((1, "E2"), (1, "E3"), (-1, "A_C")),
           "I3_gf1_plus_E1": rk((1, "D_gf1"), (1, "E1")),
           "total": rk((1, "D_EH1"), (1, "D_EH2"), (1, "D_EH3"), (1, "D_gf1"), (1, "D_gf2"), (1, "D_gf3"),
                       (1, "E1"), (1, "E2"), (1, "E3"), (-1, "A_K"), (-1, "A_C"))}
    res["seconds"] = round(time.time() - t0, 1)
    (R / "build/logs/r03_legs.json").write_text(json.dumps(res, indent=1) + "\n")
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
