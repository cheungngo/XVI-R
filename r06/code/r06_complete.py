#!/usr/bin/env python3
"""XVI-R R06 (B): Y-space completeness.
The R02-R05 Y bases (evfit.monomials) contain only EVANESCENT monomials (>= 1 hatted contraction).  The complete space
of local gh -1 three-point functionals at the power-counting dimension also contains the purely 4D monomials (only bar
contractions).  Their sigma-bar images can carry evanescent parts (d-dimensional kinetic terms, traces, ghost-eq T map).
Here: (1) enumerate the 4D monomials of every family with an independent enumerator and check that ev + 4D = all
pairings (count identity); (2) build sigma-bar of every 4D generator, split each image into evanescent coordinates and a
4D remainder (in a 4D monomial basis of each X sector); (3) the combinations with ZERO 4D remainder are additional
evanescent-trivial directions; add them and recompute the hat-aware reduced cohomology; (4) stability of the 18."""
import itertools, json, pickle, sys, time  # noqa: E401
from pathlib import Path
import sympy as sp
from sympy import QQ_I
from sympy.polys.matrices import DomainMatrix
R = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(R / "r05/code"))
import r05_final as F5  # noqa: E402
F4 = F5.F4
F = F4.F
S = F.S
I = sp.I
OUT = R / "build/logs"
DM, rk, hcat = F.DM, F.rk, F.hcat
MOM = ("p1", "p2")


def enum(pols, nmom, mode):
    """mode 'all' | '4d' | 'ev': all pairings of polarisation slots with each other / with p1,p2, H or T between
    polarisations (T never with a momentum), remaining momentum degree filled by p.p (H only)."""
    slots = [n for n in sorted(pols) for _ in range(pols[n])]
    out = set()

    def rec(rem, ml, acc, nT):
        if not rem:
            if ml % 2 == 0 and (mode == "all" or (mode == "4d" and nT == 0) or (mode == "ev" and nT > 0)):
                for pp in itertools.combinations_with_replacement([S.U.H("p1", "p1"), S.U.H("p1", "p2"), S.U.H("p2", "p2")], ml // 2):
                    out.add(acc * sp.Mul(*pp))
            return
        x, rest = rem[0], list(rem[1:])
        if ml > 0:
            for m in MOM:
                rec(rest, ml - 1, acc * S.U.H(x, m), nT)
        seen = set()
        for j, y in enumerate(rest):
            if y in seen:
                continue
            seen.add(y)
            r2 = rest[:j] + rest[j + 1:]
            rec(r2, ml, acc * S.U.H(x, y), nT)
            rec(r2, ml, acc * S.U.T(x, y), nT + 1)
    rec(slots, nmom, sp.Integer(1), 0)
    return sorted(out, key=str)


SECT = {"K": ({"aK": 2, "gc1": 1, "a3": 2}, 3), "C": ({"gb2": 1, "gc1": 1, "a3": 2}, 4),
        "S": ({"z": 1, "gc1": 1, "gc2": 1}, 3), "H": ({"a1": 2, "a2": 2, "a3": 2}, 4)}
FAMS = {"F": ({"aK": 2, "a2": 2, "a3": 2}, 2), "G": ({"gv": 1, "gc1": 1, "a3": 2}, 2), "KKc": ({"aK": 2, "a3": 2, "gc1": 1}, 1)}


def split(expr, sec, M4):
    v, rem = S.triv.coords(sp.expand(expr), S.MS[sec], sec)
    w, rem2 = S.triv.coords(sp.expand(rem), M4[sec], sec) if rem != 0 else ([0] * len(M4[sec]), 0)
    return v, w, rem2


def images_4d(fam, m):
    """sector expressions of sigma-bar(Y) for a 4D generator m (same constructions as R03/R04/R05)."""
    if fam == "F":
        f = sp.expand(m + S.R.remap(m, S.SW23))
        if f == 0:
            return None
        xh = 0
        for i in "123":
            fi = S.R.remap(f, S.REL[i]) if S.REL[i] else f
            fi = S.R.polarise(fi, "aK", "KU", "KW")
            eh, gf = S.G.S2_EH_gf(S.hhc.ext(f"a{i}"), S.G.probe_mat("E"), S.MOM[i])
            xh += S.R.strip_t(S.R.contract_id(fi, I * eh, "K", "E", "h"))
        xk = S.A_on_K(f)
        xc = F.partner_C(xk)
        return {"H": sp.expand(xh), "K": -sp.expand(xk), "C": -sp.expand(xc)}
    if fam == "G":
        xk = sp.expand(-I * S.U.H("aK", "p2") * S.hhc.subst_vec(m, "gv", {"aK": 1}))
        xc = sp.expand(S.RHO * S.Tmap(sp.nsimplify(xk, rational=True))) if xk != 0 else 0
        s = S.R.remap(m, {"p1": S.R.P1, "p2": S.R.P3, "gv": "z"})
        s = S.R.strip_t(-I * S.R.polarise(sp.expand(s), "a3", "p2", "gc2"))
        xs = sp.expand(s - S.R.remap(s, {"p1": S.R.P2, "p2": S.R.P1, "gc1": "gc2", "gc2": "gc1"}))
        return {"K": xk, "C": xc, "S": -xs}
    if fam == "KKc":
        sw = {"p1": S.R.P1, "p2": S.R.P3, "aK": "a3", "a3": "aK"}
        f = sp.expand(m - S.R.remap(m, sw))
        if f == 0:
            return None
        fp = S.R.polarise(f, "a3", "KU", "KW")
        eh, gf = S.G.S2_EH_gf(S.hhc.ext("a3"), S.G.probe_mat("E"), S.hhc.p3)
        xk = S.R.strip_t(S.R.contract_id(fp, I * eh, "K", "E", "h"))
        return {"K": xk, "C": F.partner_C(xk)}


def main():
    t0 = time.time()
    S.SIG = sp.sympify(json.loads((S.LOG / "q07_hhc_tree.json").read_text())["sigma_F"])
    import triv  # noqa: F401
    S.triv = sys.modules["triv"]
    res = {"counts": {}}
    # (1) enumeration identity, for Y families and X sectors
    for name, (p, k) in list(FAMS.items()) + [(f"X_{s}", v) for s, v in SECT.items()]:
        a, e4, ev = enum(p, k, "all"), enum(p, k, "4d"), enum(p, k, "ev")
        frozen = S.monomials(p, k)
        res["counts"][name] = {"all": len(a), "4d": len(e4), "ev": len(ev), "frozen_ev": len(frozen),
                               "frozen_equals_independent_ev": set(map(str, frozen)) == set(map(str, ev)),
                               "partition_ok": len(a) == len(e4) + len(ev)}
    M4 = {s: enum(p, k, "4d") for s, (p, k) in SECT.items()}
    off4, o = {}, 0
    for s in "KCSH":
        off4[s] = o
        o += len(M4[s])
    n4 = o
    n = S.NDIM
    # (2) images of 4D generators
    EV, FD, tags, leak = [], [], [], 0
    for fam, (p, k) in FAMS.items():
        for m in enum(p, k, "4d"):
            im = images_4d(fam, m)
            if im is None:
                continue
            ev, fd = [0] * n, [0] * n4
            for sec, ex in im.items():
                if ex == 0:
                    continue
                v, w, rem = split(ex, sec, M4)
                leak += rem != 0
                for i, c in enumerate(v):
                    ev[S.OFF[sec] + i] += c
                for i, c in enumerate(w):
                    fd[off4[sec] + i] += c
            if any(sp.expand(x) != 0 for x in ev + fd):
                EV.append(ev)
                FD.append(fd)
                tags.append(fam)
    res["n_4d_generators_nonzero"] = {f: tags.count(f) for f in FAMS}
    res["leak_outside_ev_plus_4d"] = leak
    EVm, FDm = DM(sp.Matrix(EV).T), DM(sp.Matrix(FD).T)
    Z = FDm.nullspace().transpose()                       # combos with zero 4D part
    Tev = EVm * Z if Z.shape[1] else DomainMatrix.zeros((n, 0), QQ_I)
    res["dim_4d_combos_pure_evanescent"] = Z.shape[1]
    res["rank_pure_ev_images"] = rk(Tev) if Tev.shape[1] else 0
    res["rank_ev_part_all_4d_images"] = rk(EVm)
    # (3) hat-aware reduced system
    d = pickle.loads((S.LOG / "q08_triv.pkl").read_bytes())
    B = pickle.loads((S.LOG / "q11_blocks.pkl").read_bytes())
    NB = pickle.loads((OUT / "r05_blocks.pkl").read_bytes())
    Ag = DM(d["Ag"])
    Ab = DM(pickle.loads((OUT / "r03_reduced_hhc.pkl").read_bytes()))
    Abar = Ab.vstack(DM(NB["ccK"])).vstack(DM(NB["ccb"])).vstack(DM(B["bose"]))
    AAb = Abar.vstack(Ag)
    if Tev.shape[1]:
        P = AAb * Tev
        res["pure_ev_images_inside_Nbar"] = P.is_zero_matrix
    sb = pickle.loads((OUT / "r03_sigma_bar2.pkl").read_bytes())
    cols = [v for v, t in zip(sb["cols"], sb["tags"]) if t in ("F", "KKc")]
    for m in S.monomials({"gv": 1, "gc1": 1, "a3": 2}, 2):
        v = F4.G_image(m, -1)
        if any(sp.expand(x) != 0 for x in v):
            cols.append(v)
    cl0, Nb, T0, Bg = F5.classify(AAb, Abar, Ag, cols, n)
    res["R05_reproduced"] = cl0
    T1 = hcat(T0, Tev) if Tev.shape[1] else T0
    W0, W1 = rk(hcat(T0, Bg)), rk(hcat(T1, Bg))
    res["enlarged"] = {"dim_T": rk(T1), "dim_T_plus_blind": W1, "dim_H": Nb.shape[1] - W1,
                       "new_trivial_directions": W1 - W0}
    # (4) transport of R05's 18
    lab = [(sec, mm) for sec in "KCSH" for mm in S.MS[sec]]
    idx = {(sec, F4.tex_mono(mm)): i for i, (sec, mm) in enumerate(lab)}
    reps = []
    for b in json.loads((R / "certificates/R05_H_basis.json").read_text())["basis"]:
        v = [0] * n
        for sec, mm, c in b["terms"]:
            v[idx[(sec, mm)]] = F4.coef(c)
        reps.append(v)
    C = DM(sp.Matrix(reps).T)
    Wt = hcat(T1, Bg)
    res["transport_R05_18"] = {"dim": rk(C), "trivial_after_enlargement": rk(C) + W1 - rk(hcat(C, Wt)),
                               "stable": rk(C) + W1 - rk(hcat(C, Wt)) == 0}
    res["seconds"] = round(time.time() - t0, 1)
    (OUT / "r06_complete.json").write_text(json.dumps(res, indent=1, default=str) + "\n")
    print(json.dumps(res, indent=1, default=str))


if __name__ == "__main__":
    main()
