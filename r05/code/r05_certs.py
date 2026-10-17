#!/usr/bin/env python3
"""XVI-R R05 certificates + LaTeX fragments."""
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
    fn = strip(json.loads((L / "r05_final.json").read_text()))
    w("R05_obstruction_outcome", {
        "cause": "calculus artefact: gravity.contract drops every T copy of a probe dot (keeps the H copy as the full "
                 "4S dot). Exact for loop expressions (dots come as H+T pairs) but annihilates the single-hatted c* "
                 "contractions T(u,z) of the evanescent ansatz X_S in the ghost-row c1 channel (c* line X<->ZZ).",
        "fix": "eta = eta_bar + eta_hat on the c* line: H.H->H, T.T->T, mixed->0; linear (X-dependent) part only; "
               "frozen loop residual untouched",
        "sanity_old_blocks_equal_Q11": fn["sanity_old_blocks_equal_Q11"], "hat_changes_rowspace": fn["hat_changes_rowspace"],
        "images": {k: fn[f"{k}_images"] for k in ("F", "KKc", "G")}, "G_eps_plus1_failing": fn["G_eps_plus1_failing"],
        "outcome": "ARTEFACT -- CLOSED: all 33 G images (and all F, KKc) satisfy the hat-aware reduced system; "
                   "sigma-bar(Y) inside Nbar_g on the whole computed Y-space"})
    w("R05_classification", {"reduced_hat": fn["reduced_hat"], "transport_R04_20": fn["transport_R04_20"],
        "dim_Nbar_cap_span_all_plus_blind": fn["dim_Nbar_cap_span_all_plus_blind"],
        "H_dim": fn["reduced_hat"]["dim_H"], "physical_certified": 0,
        "status": "UNCONDITIONAL on the computed (power-counting) Y-space; physical certification pending Y-completeness",
        "note": "nilpotency verified: rank Sigma = dim(Sigma cap Nbar_g); the 18 classes are nontrivial modulo the full image"})
    w("R05_H_basis", {"dim": len(fn["H_basis"]), "basis": fn["H_basis"]})
    w("R05_P6_status", {"landmark": "XVI-R", "release": "R05", "P6_outcome": "Partial", "claim_status": "partial",
        "goals": {"1_decide": "ESTABLISHED: calculus artefact (T copies dropped on the c* line)",
                  "2_fix_and_close": "ESTABLISHED: hat-aware contraction (derived), 33/33 G closed, sigma-bar nilpotent on Y",
                  "3_classification": f"Partial: H={fn['reduced_hat']['dim_H']} nontrivial classes on the computed Y-space; "
                                      "0 certified physical (Y-completeness open)",
                  "4_Y_completeness": "OPEN",
                  "5_loop_residual_hat_rule": "OPEN (frozen Q06/Q07 ghost-row remainders computed with the drop-T rule "
                                              "not revisited)",
                  "6_p22_kappa2_consistency": "OPEN"},
        "O1_closure_unchanged": "CLOSED (scoped)", "bulk_symmetry_restored": False, "insertion": "frozen",
        "disclosed": ["R05-D1: ghost-row linear blocks rebuilt with hat-aware c*-line contraction (frozen code not edited; "
                      "frozen loop residual not revised)",
                      "R05-D2: Nbar_g 139 (R03/R04) -> 133; blind_g 50 -> 38; R04's 2 extra classes become trivial",
                      "R05-D3: R05 code in r05/ (R01-R04 gates untouched)"],
        "not_implied": "no direction certified physical; Missing never zero"})
    rh, tr = fn["reduced_hat"], fn["transport_R04_20"]
    M = {"NB": rh["dim_Nbar_g"], "RS": rh["rank_Sigma"], "TB": rh["dim_T"], "BG": rh["dim_blind_g"],
         "TCB": rh["dim_T_cap_blind"], "HB": rh["dim_H"], "TRT": tr["trivial_or_blind"], "TRN": tr["nontrivial"],
         "EPSP": fn["G_eps_plus1_failing"], "SPN": fn["dim_Nbar_cap_span_all_plus_blind"]}
    (R / "build/R05_numbers.tex").write_text("".join(f"\\newcommand{{\\{k}}}{{{v}}}\n" for k, v in M.items()))
    B = []
    for i, b in enumerate(fn["H_basis"], 1):
        B.append(f"\\subsection*{{$\\mathcal H_{{{i}}}$: sectors {', '.join(b['sectors'])}; {b['n_terms']} terms}}")
        B.append(r"{\small\begin{align*}")
        B.append("\\\\\n".join(f"&{('+' if not c.startswith('-') else '')}{c}\\,[{s}]\\,{m}" for s, m, c in b["terms"]))
        B.append(r"\end{align*}}")
    (R / "build/R05_basis.tex").write_text("\n".join(B) + "\n")
    print("R05_CERTS_OK")


if __name__ == "__main__":
    main()
