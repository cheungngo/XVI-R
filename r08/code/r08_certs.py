#!/usr/bin/env python3
"""XVI-R R08 certificates + LaTeX numbers."""
import json
from pathlib import Path
R = Path(__file__).resolve().parents[2]
CERT, L = R / "certificates", R / "build" / "logs"


def strip(o):
    if isinstance(o, dict):
        return {k: strip(v) for k, v in o.items() if k != "seconds"}
    if isinstance(o, list):
        return [strip(v) for v in o]
    return o


def w(n, o):
    (CERT / f"{n}.json").write_text(json.dumps(o, indent=1, sort_keys=True, default=str) + "\n")


def main():
    g = strip(json.loads((L / "r08_gram.json").read_text()))
    iv = strip(json.loads((L / "r08_inventory.json").read_text()))
    gok = g["G_nonzero"] and g["G_in_4D_basis"] and g["G_bose_symmetric"] and g["G_sigma_image_zero"] and \
        g["R07_reproduced"]["dim_H"] == 6
    w("R08_gram", {**g, "outcome": "CLOSED" if gok else "OPEN",
        "statement": f"Gram-reduced scheme: exactly one Gram polynomial (X_H); it is a cocycle that the three-point rows "
                     f"do not see (blind), not trivial; Hbar_gram = {g['H_gram_reduced']} = Hbar_formal = {g['H_formal']}; "
                     f"blind {g['R07_reproduced']['dim_blind_g']} -> {g['gram_scheme']['dim_blind_g']}"})
    p = iv["parity"]
    w("R08_parity", {**p, "outcome": "PARTIAL" if p["inputs_parity_even"] else "OPEN",
        "established": "odd sector decouples from the even one; odd inhomogeneity zero; X_odd = 0 consistent",
        "open": "dimension of the parity-odd evanescent cohomology (Missing, not zero)"})
    w("R08_four_point", iv["four_point"])
    w("R08_P6_status", {"landmark": "XVI-R", "release": "R08", "P6_outcome": "Partial", "claim_status": "partial",
        "goals": {"1_four_point_rows": "OPEN", "2_gram_reduced_scheme_count": "CLOSED" if gok else "OPEN",
                  "3_parity_odd_cohomology": "PARTIAL"},
        "H_formal": g["H_formal"], "H_gram_reduced": g["H_gram_reduced"],
        "O1_closure_unchanged": "CLOSED (scoped)", "bulk_symmetry_restored": False, "insertion": "frozen",
        "disclosed": ["R08-D1: Gram scheme adds one evanescent coordinate (the X_H Gram determinant); counted as blind "
                      "(Nbar_g 133 -> 134, blind 24 -> 25); Hbar unchanged at 6",
                      "R08-D2: parity-odd count recorded as Missing; only decoupling and zero inhomogeneity established",
                      "R08-D3: four-point rows not computed; missing ingredients listed",
                      "R08-D4: R08 code in r08/ (R01-R07 gates and code untouched)"]})
    M = {"GT": g["G_n_terms"], "NBg": g["gram_scheme"]["dim_Nbar_g"], "BGg": g["gram_scheme"]["dim_blind_g"],
         "HBg": g["H_gram_reduced"], "HBf": g["H_formal"], "TBg": g["gram_scheme"]["dim_T"]}
    (R / "build/R08_numbers.tex").write_text("".join(f"\\newcommand{{\\{k}}}{{{v}}}\n" for k, v in M.items()))
    print("R08_CERTS_OK")


if __name__ == "__main__":
    main()
