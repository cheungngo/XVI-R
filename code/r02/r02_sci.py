#!/usr/bin/env python3
"""XVI-R R02: computed trivial span and cohomology quotient on N_g (kappa^3, three-point, pole level).
Y-space (gh -1, local, the full monomial bases used since Q08): F = Ktilde.F(h,h) (momentum degree 2, Ktilde = K + cbar
F'(h) so the cbar partner is rho*T), G = c*.G(c,h) (degree 2), P = cbar.Phi(h,h) (degree 3, standalone cbar; the ghost
equation is then imposed on the image, not assumed).  sigma is the linear map Y -> X built with the Q11 channel operators
and derived signs (anti - field).  Computed (not hand-picked):
  Z   = ker([A;Ag] Sigma),  T = Sigma Z  (sigma-images that satisfy every constraint: the trivial span in N_g),
  W'  = T + blind_g,  H_cand = N_g / W'.
Robust test: a quotient direction outside span(Sigma) + blind_g (ignoring constraints) cannot be reached by ANY combination
of the computed images, even failing ones.  Nothing here proves nontriviality in cohomology: the Q11 composition gap is
open (sigma(Y) is not inside N), so every non-trivial class is a CANDIDATE."""
import json, pickle, sys, time  # noqa: E401
from pathlib import Path
import sympy as sp
from sympy import QQ_I
from sympy.polys.matrices import DomainMatrix
R = Path(__file__).resolve().parents[2]
QQ = R / "build" / "xviq" / "XVI-Q"
sys.path.insert(0, str(QQ / "code/p6q06/q07/q08/q09/q10/q11"))
import sigma11 as S  # noqa: E402
OUT = R / "build" / "logs"


def DM(M):
    return DomainMatrix.from_Matrix(M).convert_to(QQ_I)


def rk(D):
    return D.rank() if D.shape[0] and D.shape[1] else 0


def hcat(*Ds):
    Ds = [d for d in Ds if d.shape[1]]
    return Ds[0].hstack(*Ds[1:]) if len(Ds) > 1 else Ds[0]


def colbasis(D):
    """independent columns of D (pivot columns)."""
    if not D.shape[1]:
        return D
    _, piv = D.transpose().transpose().rref()
    return D.extract(list(range(D.shape[0])), list(piv))


def family_P():
    out = []
    for m in S.monomials({"gv": 1, "a2": 2, "a3": 2}, 3):
        f = sp.expand(m + S.R.remap(m, S.SW23))
        if f == 0:
            continue
        vh, l1 = S.vec("H", S.E_op(f))
        vc, l2 = S.vec("C", S.A_on_C(f))
        out.append({"anti": vh, "field": vc, "leak": l1 or l2})
    return out


def tex_mono(m):
    parts = []
    for b, e in sp.sympify(m).as_powers_dict().items():
        t, x, y = b.name.split("_")
        nm = {"aK": "K", "gc1": "c_1", "gc2": "c_2", "gb2": r"\bar c", "z": "c^*"}

        def f(s):
            if s in nm:
                return nm[s]
            return (r"\epsilon_" if s[0] == "a" else "p_") + s[1:]
        s = (r"(\hat{%s}\!\cdot\!\hat{%s})" if t == "T" else r"(%s\!\cdot\!%s)") % (f(x), f(y))
        parts.append(s + (f"^{e}" if e > 1 else ""))
    return "".join(parts)


def main():
    t0 = time.time()
    S.SIG = sp.sympify(json.loads((S.LOG / "q07_hhc_tree.json").read_text())["sigma_F"])
    d = pickle.loads((S.LOG / "q08_triv.pkl").read_bytes())
    A, Ag = d["A"], d["Ag"]
    n = A.shape[1]
    fams = {"F": S.family_F(), "G": S.family_G(), "P": family_P()}
    cols, tags = [], []
    for k, fam in fams.items():
        for g in fam:
            v = S.add(g["anti"], S.scale(-1, g["field"]))
            if any(sp.expand(x) != 0 for x in v):
                cols.append(v)
                tags.append(k)
    leak = {k: sum(g["leak"] for g in f) for k, f in fams.items()}
    Sig = DM(sp.Matrix(cols).T)                                           # n x nY
    AA = DM(A.col_join(Ag))
    rS = rk(Sig)
    Z = (AA * Sig).nullspace().transpose()                                # nY x z
    T = colbasis(Sig * Z) if Z.shape[1] else DomainMatrix.zeros((n, 0), QQ_I)
    ZA = (DM(A) * Sig).nullspace().transpose()
    TA = Sig * ZA if ZA.shape[1] else DomainMatrix.zeros((n, 0), QQ_I)
    by_fam = {}
    for k in fams:
        idx = [j for j, t in enumerate(tags) if t == k]
        Sk = Sig.extract(list(range(n)), idx)
        Zk = (AA * Sk).nullspace().transpose()
        by_fam[k] = {"n_generators": len(idx), "rank_sigma": rk(Sk), "dim_trivial_in_Ng": rk(Sk * Zk) if Zk.shape[1] else 0}
    Ng = AA.nullspace().transpose()                                       # n x 104
    zc = [j for j in range(n) if all(A[i, j] == 0 for i in range(A.shape[0]))]
    Bm = DM(sp.Matrix([[1 if i == j else 0 for j in zc] for i in range(n)]))
    kb = (DM(Ag) * Bm).nullspace().transpose()
    Bg = Bm * kb
    # Q11 W (individually passing images) for the 43 baseline
    keep = pickle.loads((S.LOG / "q11_sigma_vecs.pkl").read_bytes())
    W11 = hcat(DM(sp.Matrix([*keep[("F", "derived_minus1")], *keep[("G", "derived_minus1")]]).T), Bg)
    dNg, dW11 = Ng.shape[1], rk(W11)
    Wp = hcat(W11, T, Bg)
    dWp = rk(Wp)
    Wrob = hcat(Sig, Bg)
    # adapted basis of the 43: N_g directions reduced modulo W11 -> split into (in W'), (in span Sig+blind), rest
    res = {"n_coords": n, "leak": leak, "n_Y_generators": len(cols), "rank_Sigma": rS,
           "dim_kernel_constraints_on_Y": Z.shape[1], "dim_T_trivial_in_Ng": rk(T),
           "dim_T_rows_only": rk(TA) if TA.shape[1] else 0, "by_family": by_fam,
           "dim_N_g": dNg, "dim_W_Q11": dW11, "unclassified_Q11": dNg - dW11,
           "dim_Wprime": rk(hcat(Wp, Ng)) and dWp, "dim_Wprime_cap_Ng": rk(hcat(Wp)) + dNg - rk(hcat(Wp, Ng)),
           "dim_Ng_plus_Wprime": rk(hcat(Wp, Ng)),
           "dim_Ng_cap_spanSigma_plus_blind": rk(Wrob) + dNg - rk(hcat(Wrob, Ng))}
    # quotient basis: greedy over N_g basis columns, add those independent of current span (W' cap N_g handled by rank)
    cur = hcat(Wp)
    base = rk(hcat(Wp, Ng)) - dNg                                          # dims of Wp outside N_g
    reps = []
    for j in range(dNg):
        c = Ng.extract(list(range(n)), [j])
        t = hcat(cur, c)
        if rk(t) > rk(cur):
            cur = t
            reps.append(c)
    res["dim_H_candidate"] = len(reps)
    # robust: which reps remain outside span(Sig)+blind (sequentially)
    curR = hcat(Wrob, W11)
    rob = 0
    for c in reps:
        t = hcat(curR, c)
        if rk(t) > rk(curR):
            curR = t
            rob += 1
    res["dim_H_candidate_robust"] = rob
    res["base_note"] = base
    # 43 classification counts (adapted basis)
    newly = (dNg - dW11) - len(reps)
    res["classification_43"] = {"trivial_computed_sigma": newly, "blind": 0,
                                "cohomology_candidate": len(reps), "of_which_robust": rob, "physical_proven": 0}
    # canonical representatives: rref of W' (pivots P); reduce N_g modulo W' to vanish on P; rref the residual.
    # Unique for the fixed coordinate order (sectors K, C, S, H; monomial order of Q08), hence basis-independent.
    Wr, Pw = Wp.transpose().to_field().rref()
    Wr = Wr.to_Matrix()[:len(Pw), :]
    Nr = Ng.transpose().to_field().to_Matrix()
    for i, pc in enumerate(Pw):
        for r in range(Nr.shape[0]):
            if Nr[r, pc] != 0:
                Nr[r, :] = Nr[r, :] - Nr[r, pc] * Wr[i, :]
    Hr, Ph = DM(Nr).to_field().rref()
    Hr = Hr.to_Matrix()[:len(Ph), :]
    assert len(Ph) == len(reps)
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
    res["basis"] = basis
    res["seconds"] = round(time.time() - t0, 1)
    (OUT / "r02_sci.json").write_text(json.dumps(res, indent=1, default=str) + "\n")
    print(json.dumps({k: v for k, v in res.items() if k != "basis"}, indent=1, default=str))


if __name__ == "__main__":
    main()
