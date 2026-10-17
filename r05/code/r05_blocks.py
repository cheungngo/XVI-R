#!/usr/bin/env python3
"""XVI-R R05: hat-aware contraction over the c* line in the ghost-row c1 channel.
gravity.contract (frozen) replaces a 4S dot by its H copy and DROPS every T copy of a probe dot.  This is exact for loop
expressions whose probe dots always come as H+T pairs, but not for the evanescent ansatz X_S = Gamma_{c*cc} (pole CT),
whose monomials carry a single hatted c* contraction T(u, z) with no H partner: those terms are annihilated.
With eta = eta_bar + eta_hat the correct rule for the c* line (probes X in the kinetic, ZZ in X_S) is
    H_X_a H_ZZ_b -> H_a_b,  T_X_a T_ZZ_b -> T_a_b,  mixed -> 0 .
Only the LINEAR (X-dependent) part of the ccK/ccb rows is rebuilt; the frozen loop residual (base) is not touched,
so N-type kernels (homogeneous) are affected and no frozen certificate is revised."""
import pickle, sys, time  # noqa: E401
from pathlib import Path
import sympy as sp
R = Path(__file__).resolve().parents[2]
Q = R / "build/xviq/XVI-Q/code"
sys.path[:0] = [str(Q / "p6q05"), str(Q / "p6q06"), str(Q / "p6q06/q07")]
import gravity as G  # noqa: E402
import rows as RW  # noqa: E402
import evfit2 as E2  # noqa: E402
from evfit import ansatz  # noqa: E402
OUT = R / "build/logs"
_orig = G.contract


def contract_hat(expr, pairs):
    if [tuple(p) for p in pairs] != [("X", "ZZ")]:
        return _orig(expr, pairs)
    expr = sp.expand(expr)
    out = 0
    for term in sp.Add.make_args(expr):
        xs, zs, coef = None, None, term
        for s, e in term.as_powers_dict().items():
            if not isinstance(s, sp.Symbol):
                continue
            p = s.name.split("_")
            if len(p) != 3:
                continue
            if "X" in p[1:]:
                xs = (p[0], p[2] if p[1] == "X" else p[1], e)
                coef = coef / s ** e
            elif "ZZ" in p[1:]:
                zs = (p[0], p[2] if p[1] == "ZZ" else p[1], e)
                coef = coef / s ** e
        if xs is None or zs is None:
            out += _orig(term, pairs)
            continue
        assert xs[2] == 1 and zs[2] == 1
        if xs[0] != zs[0]:
            continue
        a, b = sorted([xs[1], zs[1]])
        out += coef * sp.Symbol(f"{xs[0]}_{a}_{b}")
    return sp.expand(out)


def lin_block(target, hat):
    XK, xk = ansatz({"aK": 2, "gc1": 1, "a3": 2}, 3, "xk")
    XC, xc = ansatz({"gb2": 1, "gc1": 1, "a3": 2}, 4, "xc")
    XS, xs = ansatz({"z": 1, "gc1": 1, "gc2": 1}, 3, "xs")
    XH, xh = ansatz({"a1": 2, "a2": 2, "a3": 2}, 4, "xh")
    xall = xk + xc + xs + xh
    G.contract = contract_hat if hat else _orig
    try:
        XA = XK if target == "ccK" else XC
        r = E2.ghost_row_total(target, XA, XS) - E2.ghost_row_total(target, sp.Integer(0), sp.Integer(0))
    finally:
        G.contract = _orig
    r = sp.expand(r)
    X = set(xall)
    gens = sorted([s for s in r.free_symbols if s not in X], key=str)
    return sp.linear_eq_to_matrix(sp.Poly(r, *gens).coeffs(), xall)[0]


def canon(sym):
    return sym


def main():
    t0 = time.time()
    out = {}
    for t in ("ccK", "ccb"):
        out[t] = lin_block(t, True)
        out[t + "_old"] = lin_block(t, False)
    pickle.dump(out, open(OUT / "r05_blocks.pkl", "wb"))
    print({k: v.shape for k, v in out.items()}, round(time.time() - t0, 1))


if __name__ == "__main__":
    main()
