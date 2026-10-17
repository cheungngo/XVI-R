#!/usr/bin/env python3
"""XVI-R R07: hat-aware contraction on the EXTERNAL graviton-probe line (K~ <-> E) used in row assembly
(rows.contract_id(..., 'K', 'E', 'h')), the analogue of R05's c*-line rule.  The frozen gravity.contract drops every
T copy of a probe dot and keeps the H copy as the full 4S dot; exact only if BOTH sides carry 4S dots (H+T pairs).
A K-side carrying a hat-only (ghat) or tilde-only (gt) contraction -- evanescent ansatz monomials T(aK, .), the D13
E2 structure, ghat from symmetric integration -- is mistreated.  Rule (eta = eta_bar + eta_hat): H.H -> H, T.T -> T,
mixed -> 0.  Applied ONLY to pairs whose probes are all external-line probes {KU,KW,EU,EW}; internal loop lines untouched.
This module (i) measures the effect on the reduced hhc linear block and on the loop channels C and D_EH."""
import json, pickle, sys, time  # noqa: E401
from pathlib import Path
import sympy as sp
R = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(R / "r06/code"))
sys.path.insert(0, str(R / "r07/code"))
sys.path.insert(0, str(R / "r03/code"))
import r06_refit as F6  # noqa: E402
B5, G, S = F6.B5, F6.G, F6.S
EXTP = {"KU", "KW", "EU", "EW", "X", "ZZ"}
OUT = R / "build/logs"


def _pd(s):
    p = s.name.split("_")
    return (p[0], p[1], p[2]) if len(p) == 3 and p[0] in ("H", "T") else None


def contract_hatgen(expr, pairs):
    pairs = [tuple(p) for p in pairs]
    if not all(a in EXTP and b in EXTP for a, b in pairs):
        return B5._orig(expr, pairs)
    expr = sp.expand(expr)
    out = 0
    for term in sp.Add.make_args(expr):
        coef, ok, dots = term, True, {}
        for s, e in term.as_powers_dict().items():
            if not isinstance(s, sp.Symbol):
                continue
            pd = _pd(s)
            if pd is None:
                if s.name.startswith("s_") and any(x in s.name.split("_")[1:] for pr in pairs for x in pr):
                    ok = False
                continue
            k, a, b = pd
            allp = {x for pr in pairs for x in pr}
            inv = [x for x in (a, b) if x in allp]
            if not inv:
                continue
            if len(inv) > 1 or e != 1 or inv[0] in dots:
                ok = False
                continue
            pn = inv[0]
            dots[pn] = (k, b if a == pn else a)
            coef = coef / s
        if not ok or any(p not in dots or q not in dots for p, q in pairs):
            out += B5._orig(term, pairs)
            continue
        val = coef
        for p, q in pairs:
            (k1, x), (k2, y) = dots[p], dots[q]
            if k1 != k2:
                val = 0
                break
            if k1 == "T" and (x in ("p1", "p2") or y in ("p1", "p2")):
                val = 0
                break
            val = val * (sp.Symbol(f"{k1}_{min(x, y)}_{max(x, y)}") if x != y else sp.Symbol(f"{k1}_{x}_{x}"))
        out += val
    return sp.expand(out)


def with_hat(fn, *a):
    G.contract = contract_hatgen
    try:
        return fn(*a)
    finally:
        G.contract = B5._orig


def main():
    t0 = time.time()
    import r03_fix as F3
    res = {}
    old = pickle.loads((OUT / "r03_reduced_hhc.pkl").read_bytes())
    new = with_hat(F3.reduced_hhc_block)
    pickle.dump(new, open(OUT / "r07_reduced_hhc_hat.pkl", "wb"))
    o, n = F6.rk(old), F6.rk(new)
    res["reduced_hhc_block"] = {"rank_old": o, "rank_new": n, "rank_joint": F6.rk(old.col_join(new)),
                                "rowspace_changed": F6.rk(old.col_join(new)) != o or o != n}
    # loop channels with external-probe contraction
    import r07_hhc as H7
    hhc = S.hhc
    cO, cN = hhc.chan_C(), with_hat(hhc.chan_C)
    res["chan_C_changed"] = sp.expand(cO - cN) != 0
    res["chan_C_diff_4D"] = str(F6.t0proj(sp.expand(cN - cO)))
    dd = 0
    for i in ("2", "3"):
        dO = H7.chan_D_split(i)[0]
        dN = with_hat(H7.chan_D_split, i)[0]
        dd += sp.expand(dN - dO)
    res["chan_D_EH_changed"] = sp.expand(dd) != 0
    res["chan_D_EH_diff_4D"] = str(F6.t0proj(dd))
    res["seconds"] = round(time.time() - t0, 1)
    (OUT / "r07_hatgen.json").write_text(json.dumps(res, indent=1) + "\n")
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
