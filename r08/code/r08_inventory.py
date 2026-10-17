#!/usr/bin/env python3
"""XVI-R R08 fronts 1 and 3 (honest inventory + decoupling proof pieces).
Front 3 (parity-odd): every input of the kappa^3 system (stored loop poles, tree vertices, ansatz monomials, the R07
corrected row blocks' generating expressions, p22 two-point poles) is built from H_, T_ (and loop s_) dots only; the
calculus has no Levi-Civita object.  sigma-bar acts by linear substitutions of vectors and metric contractions, which
preserve the number of epsilon tensors mod 2.  Hence the parity-odd sector decouples, its inhomogeneity is ZERO, and
X_odd = 0 is a consistent choice.  The DIMENSION of the odd evanescent cohomology needs an epsilon calculus in DRED
(4S epsilon with tilde components) that does not exist in the frozen code: Missing (not zero).
Front 1 (four-point rows): records which ingredients exist in the frozen calculus."""
import json, re, sys  # noqa: E401
from pathlib import Path
R = Path(__file__).resolve().parents[2]
Q = R / "build/xviq/XVI-Q"
OUT = R / "build/logs"
ALLOWED = re.compile(r"^(H|T)_[A-Za-z0-9]+_[A-Za-z0-9]+$|^s_[A-Za-z0-9]+$|^(t|L|I|sigma_F|d|eps|U)$")


def scan(expr_str):
    import sympy as sp
    e = sp.sympify(expr_str)
    bad = sorted({s.name for s in e.free_symbols if not ALLOWED.match(s.name)})
    return bad


def main():
    res = {"parity": {}, "four_point": {}}
    files = sorted(p for pat in ("q05_*.json", "q06_*.json", "q07_*.json") for p in (Q / "build/logs").glob(pat))
    bad = {}
    for f in files:
        d = json.loads(f.read_text())
        vals = []

        def walk(o):
            if isinstance(o, dict):
                for v in o.values():
                    walk(v)
            elif isinstance(o, str) and ("H_" in o or "T_" in o):
                vals.append(o)
        walk(d)
        b = sorted({x for v in vals for x in scan(v)})
        if b:
            bad[f.name] = b
    for f in ("r07_Rbar_hatgen.txt", "r07_Xbar_hatgen.json", "r07_p22_completion.json"):
        p = OUT / f
        txt = p.read_text()
        vals = list(json.loads(txt).values()) if f.endswith(".json") else [txt]
        b = sorted({x for v in vals if v.strip() not in ("", "0") for x in scan(v)})
        if b:
            bad[f] = b
    src = "\n".join(p.read_text() for p in (Q / "code").rglob("*.py"))
    res["parity"]["non_dot_symbols_in_inputs"] = bad
    res["parity"]["levi_civita_in_frozen_calculus"] = bool(re.search(r"levi|epsilon_tensor|LeviCivita|eps4", src))
    res["parity"]["inputs_parity_even"] = not bad and not res["parity"]["levi_civita_in_frozen_calculus"]
    res["parity"]["odd_inhomogeneity_zero"] = res["parity"]["inputs_parity_even"]
    res["parity"]["odd_cohomology_dimension"] = "Missing (no DRED epsilon calculus)"
    g = (Q / "code/p6q05/gravity.py").read_text()
    res["four_point"] = {
        "V4_graviton_vertex": bool(re.search(r"def V4", g)),
        "K_c_hh_vertex": bool(re.search(r"def V_K_c_hh", g)),
        "cbar_c_hh_vertex": bool(re.search(r"def V_cb_c_hh", g)),
        "four_point_loop_pole_runs": any((Q / "build/logs").glob("*4pt*")) or any((Q / "build/logs").glob("*hhhc*")),
        "missing": ["kappa^2 four-graviton vertex", "K~ c h h and cbar c h h vertices", "one-loop four-point UV poles "
                    "(boxes, triangles with V4, bubbles) for rows c.h.h.h, c.c.cbar.h, c.c.K.h", "four-point evanescent "
                    "ansatz and its sigma-bar map"]}
    res["four_point"]["outcome"] = "OPEN (ingredients Missing; nothing assumed zero)"
    (OUT / "r08_inventory.json").write_text(json.dumps(res, indent=1) + "\n")
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
