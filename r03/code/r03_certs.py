#!/usr/bin/env python3
"""XVI-R R03 certificates + LaTeX fragments (build/R03_*.tex)."""
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


def j(n):
    return strip(json.loads((L / f"{n}.json").read_text()))


def w(n, o):
    (CERT / f"{n}.json").write_text(json.dumps(o, indent=1, sort_keys=True, default=str) + "\n")


def main():
    dg, lg, fx, f2, fn, gs = (j(x) for x in ("r03_diag", "r03_legs", "r03_fix", "r03_fix2", "r03_final", "r03_gscan"))
    w("R03_gap_diagnosis", {"channel_coefficient_scan": {k: dg[k] for k in dg if k.startswith("channel")},
        "D_split_exact": dg["D_split_exact"], "leg_identities_hhc_block": lg,
        "reading": "defect (rank 40) = A-channel on (D_gf + E) at the ghost leg: the antighost / gauge-fixing sector; "
                   "EH gauge invariance (I1) and leg commutation (I2, E vs A_C) hold. No sign/normalisation fixes it."})
    fb = fn["H_bar_basis"]
    w("R03_reduced", {"equivalence_full_vs_reduced_on_kerAg_rank": fx["equivalence_rank_full_minus_reduced_on_kerAg"],
        "lambda_test": {k: f2[k] for k in ("lambda_candidate", "sigma_F", "lambda_closes_on_kerAg")},
        "dim_N_g": fx["dim_N_g"], "dim_Nbar_g": fx["dim_Nbar_g"], "dim_N_g_cap_Nbar_g": f2["reduced"]["dim_N_g_cap_Nbar_g"],
        "F_images_failing_reduced": fx["images_failing_reduced"]["F"], "n_F": fn["n_F"], "KKc": f2["KKc"],
        "G": {"eps_minus1_failing": fn["G_eps_-1_failing"], "eps_plus1_failing": fn["G_eps_+1_failing"],
              "fail_blocks": f2["sigma_bar_fail_blocks"], "coefficient_scan_nullspace": gs["nullspace"]},
        "gap_outcome": "hhc gap CLOSED in reduced (Gamma-bar) variables: all F and KKc images satisfy the reduced system; "
                       "residual G-sector (c*.c.h) defect in ghost rows ccK/ccb remains OPEN (no coefficient fix)"})
    w("R03_classification", {"reduced": fn["reduced"], "transport_R02_22": fn["transport_22"],
        "H_bar_dim_nonblind": fn["reduced"]["dim_H_bar_nonblind"], "physical_certified": 0,
        "status": "CONDITIONAL on the G-sector (ccK/ccb) defect; robust: the classes lie outside span(all computed images) + blind"})
    w("R03_H_bar_basis", {"dim": len(fb), "basis": fb})
    w("R03_P6_status", {"landmark": "XVI-R", "release": "R03", "P6_outcome": "Partial", "claim_status": "partial",
        "goals": {"1_diagnosis": "ESTABLISHED: B-integrated-out (on-shell) BRST; defect = antighost/gauge-fixing sector at the ghost leg",
                  "2_gap": "Partial: hhc gap CLOSED in reduced variables (F, KKc); G-sector ghost-row defect OPEN",
                  "3_classification": f"Partial: of R02's 22, {fn['transport_22']['fail_reduced_constraints']} not reduced cocycles, "
                                      f"{fn['transport_22']['trivial_or_blind_in_reduced']} trivial/blind, "
                                      f"{fn['transport_22']['nontrivial_classes_in_H_bar']} conditional classes; 0 physical",
                  "4_p22_kappa2_consistency": "OPEN"},
        "O1_closure_unchanged": "CLOSED (scoped)", "bulk_symmetry_restored": False, "insertion": "frozen",
        "disclosed": ["R03-D1: N_g (Q08 full rows, 104) differs from Nbar_g (reduced, 139); intersection 101",
                      "R03-D2: R03 code lives in r03/ so the R01/R02 gates are untouched (no frozen-layer edit)"],
        "not_implied": "no direction certified physical; Missing never zero"})
    M = {"DR": 40, "EQ": fx["equivalence_rank_full_minus_reduced_on_kerAg"], "NB": fx["dim_Nbar_g"], "NI": f2["reduced"]["dim_N_g_cap_Nbar_g"],
         "NF": fn["n_F"], "NK": f2["KKc"]["blocks"].get("pass", 0), "NKG": f2["KKc"]["n"], "GM": fn["G_eps_-1_failing"], "GP": fn["G_eps_+1_failing"],
         "TB": fn["reduced"]["dim_T_bar"], "BG": fn["reduced"]["dim_blind_g"], "TCB": fn["reduced"]["dim_T_cap_blind"],
         "HB": fn["reduced"]["dim_H_bar_nonblind"], "RS": fn["reduced"]["rank_Sigma_bar"],
         "TA": fn["transport_22"]["fail_reduced_constraints"], "TT": fn["transport_22"]["trivial_or_blind_in_reduced"],
         "TN": fn["transport_22"]["nontrivial_classes_in_H_bar"], "SR": fn["reduced"]["dim_span_Sigma_plus_blind_cap_Nbar"]}
    (R / "build/R03_numbers.tex").write_text("".join(f"\\newcommand{{\\{k}}}{{{v}}}\n" for k, v in M.items()))
    legs = [("$A\\cdot D_{EH}^{(1)}$ (EH gauge invariance)", lg["I1_EH_leg1"]),
            ("$A\\cdot(D^{(2)}+D^{(3)})-D\\cdot A_K$", lg["I2_legs23_vs_A_K"]),
            ("same, EH only", lg["I2_EH_only_legs23_vs_A_K"]),
            ("$E^{(2)}+E^{(3)}-A_C$", lg["I3_E_legs23_vs_A_C"]),
            ("$A\\cdot(D_{gf}^{(1)}+E^{(1)})$", lg["I3_gf1_plus_E1"]), ("total", lg["total"])]
    (R / "build/R03_legs.tex").write_text("\\newcommand{\\LegRows}{" + "\n".join(f"{a} & {b}\\\\" for a, b in legs) + "}\n")
    B = []
    for i, b in enumerate(fb, 1):
        B.append(f"\\subsection*{{$\\bar{{\\mathcal O}}_{{{i}}}$: sectors {', '.join(b['sectors'])}; {b['n_terms']} terms}}")
        B.append(r"{\small\begin{align*}")
        B.append("\\\\\n".join(f"&{('+' if not c.startswith('-') else '')}{c}\\,[{s}]\\,{m}" for s, m, c in b["terms"]))
        B.append(r"\end{align*}}")
    (R / "build/R03_basis.tex").write_text("\n".join(B) + "\n")
    print("R03_CERTS_OK")


if __name__ == "__main__":
    main()
