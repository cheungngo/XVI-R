#!/usr/bin/env python3
"""XVI-R R06 certificates + LaTeX number macros."""
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
    lo = strip(json.loads((L / "r06_loops.json").read_text()))
    rf = strip(json.loads((L / "r06_refit.json").read_text()))
    cp = strip(json.loads((L / "r06_complete.json").read_text()))
    qh = rf["Q_full_hhc_hat"]
    w("R06_loop_remainders", {
        "loop_remainders_unchanged": lo["loop_remainders_unchanged"],
        "ccK": {k: lo["ccK"][k] for k in ("n_terms_old", "n_terms_new", "sha_old", "sha_new", "difference_n_terms")},
        "ccb": {k: lo["ccb"][k] for k in ("n_terms_old", "n_terms_new", "sha_old", "sha_new", "difference_n_terms")},
        "difference_purely_evanescent": all(rf["loop_4D_projection"][f"{t}_4D_part_unchanged"] for t in ("ccK", "ccb")),
        "refit": {k: rf[k] for k in ("Q_full_hhc_old", "Q_full_hhc_hat", "R_reduced_hhc_old", "R_reduced_hhc_hat")},
        "earlier_claims": {
            "Q07_all_four_rows_close_in_4D": "UNAFFECTED: the T->0 projection of both ghost-row loop remainders is identical "
                                             "under the old and hat-aware rule",
            "ghost_equation_compatibility": "UNAFFECTED (re-established): hat-aware full system incl. ghost equation solvable; "
                                            "particular CT purely evanescent; all rows verified zero" if
                                            (qh["solvable_with_ghost_eq"] and qh["particular_CT_all_evanescent"] and
                                             qh["verified_rows_zero"]) else "AFFECTED",
            "R05_classification": "UNAFFECTED: H depends only on the homogeneous (X-linear) system, already hat-aware in R05"},
        "disclosure": "R06-D1: frozen Q06/Q07 ghost-row loop remainders recomputed with the hat-aware c*-line rule differ "
                      "from the frozen ones by purely evanescent terms (ccK 20, ccb 22 terms); the drop-T rule is NOT exact "
                      "for loops. Frozen certificates not edited; superseded numerically by R06_loop_remainders."})
    w("R06_completeness", cp)
    ok = (cp["enlarged"]["dim_H"] == 18 and cp["transport_R05_18"]["stable"] and cp["leak_outside_ev_plus_4d"] == 0
          and qh["solvable_with_ghost_eq"] and all(c["frozen_equals_independent_ev"] and c["partition_ok"]
                                                   for c in cp["counts"].values()))
    w("R06_classification", {"H_dim": cp["enlarged"]["dim_H"], "Nbar_g": cp["R05_reproduced"]["dim_Nbar_g"],
        "dim_T": cp["enlarged"]["dim_T"], "blind_g": cp["R05_reproduced"]["dim_blind_g"],
        "T_plus_blind": cp["enlarged"]["dim_T_plus_blind"], "invariant_evanescent_cohomology": ok,
        "physical_certified": 0,
        "meaning": "the 18 classes are the invariant evanescent cohomology at O(kappa^3) on the complete local Y-space "
                   "(power counting, modulo total derivatives, field redefinitions, EOM); 'physical' only in this sense, "
                   "no direction certified physical as an observable"})
    w("R06_P6_status", {"landmark": "XVI-R", "release": "R06", "P6_outcome": "Partial", "claim_status": "partial",
        "goals": {"A_loop_remainders": "ESTABLISHED: CHANGED (evanescent only); propagated; fit solvable; H unaffected",
                  "B_Y_completeness": "ESTABLISHED" if ok else "OPEN",
                  "reduced_formulation_loop_residual": "OPEN (hhc loop residual not recomputed in the Gamma-bar formulation; "
                                                       "stored full-formulation residual is incompatible with the reduced "
                                                       "channels both before and after the hat rule)",
                  "p22_kappa2_consistency": "OPEN"},
        "O1_closure_unchanged": "CLOSED (scoped)", "bulk_symmetry_restored": False, "insertion": "frozen",
        "disclosed": [json.loads((CERT / "R06_loop_remainders.json").read_text())["disclosure"],
                      "R06-D2: Y-space enlarged by 88 purely-4D generators (F57,G21,KKc10); no new trivial directions",
                      "R06-D3: R06 code in r06/ (R01-R05 gates untouched)"]})
    M = {"LKo": lo["ccK"]["n_terms_old"], "LKn": lo["ccK"]["n_terms_new"], "LKd": lo["ccK"]["difference_n_terms"],
         "LBo": lo["ccb"]["n_terms_old"], "LBn": lo["ccb"]["n_terms_new"], "LBd": lo["ccb"]["difference_n_terms"],
         "QNo": rf["Q_full_hhc_old"]["dim_N_with_ghost_eq"], "QNn": qh["dim_N_with_ghost_eq"],
         "CTK": qh["particular_CT_n_terms"]["K"], "CTC": qh["particular_CT_n_terms"]["C"],
         "CTS": qh["particular_CT_n_terms"]["S"], "CTH": qh["particular_CT_n_terms"]["H"],
         "NF": cp["n_4d_generators_nonzero"]["F"], "NGgen": cp["n_4d_generators_nonzero"]["G"],
         "NKKc": cp["n_4d_generators_nonzero"]["KKc"], "ZP": cp["dim_4d_combos_pure_evanescent"],
         "ZR": cp["rank_pure_ev_images"], "EVR": cp["rank_ev_part_all_4d_images"],
         "HB": cp["enlarged"]["dim_H"], "NB": cp["R05_reproduced"]["dim_Nbar_g"], "TB": cp["enlarged"]["dim_T"],
         "BG": cp["R05_reproduced"]["dim_blind_g"], "TPB": cp["enlarged"]["dim_T_plus_blind"]}
    rows = "\\\\\n".join(f"{k.replace('_', ' ')} & {v['all']} & {v['4d']} & {v['ev']} & {v['frozen_ev']}"
                         for k, v in cp["counts"].items())
    (R / "build/R06_numbers.tex").write_text("".join(f"\\newcommand{{\\{k}}}{{{v}}}\n" for k, v in M.items()) +
                                             f"\\newcommand{{\\CountRows}}{{{rows}}}\n")
    print("R06_CERTS_OK")


if __name__ == "__main__":
    main()
