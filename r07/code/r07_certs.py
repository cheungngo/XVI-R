#!/usr/bin/env python3
"""XVI-R R07 certificates + LaTeX numbers."""
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


def j(n):
    return strip(json.loads((L / n).read_text()))


def main():
    h1, hg, rc, pz = j("r07_hhc.json"), j("r07_hatgen.json"), j("r07_reclass.json"), j("r07_p22.json")
    f1 = rc["front1"]
    hat = rc["hatgen"]
    w("R07_contraction_Kline", {
        "finding": "gravity.contract also mistreats the external K~<->E graviton-probe line (rows.contract_id 'K','E'): "
                   "it drops T copies and keeps H as the full 4S dot, exact only when both sides carry 4S dots. K-side "
                   "ghat/gt contractions (evanescent X_K monomials T(aK,.), D13 E2 structure) are lost or mis-hatted.",
        "rule": "eta = eta_bar + eta_hat on every external probe line: H.H->H, T.T->T, mixed->0 (R05 rule, extended)",
        "reduced_hhc_block": hg["reduced_hhc_block"], "loop_chan_C_changed": hg["chan_C_changed"],
        "loop_chan_C_change_4D": hg["chan_C_diff_4D"], "loop_chan_D_EH_changed": hg["chan_D_EH_changed"],
        "loop_chan_D_EH_change_4D": hg["chan_D_EH_diff_4D"], "ghost_blocks_equal_R05": rc["ghost_blocks_equal_R05"]})
    w("R07_classification", {"R05rule_sanity": rc["R05rule"], "hatgen": hat, "H_dim": hat["dim_H"],
        "H_basis_dim": len(rc["hatgen_H_basis_vectors"]), "physical_certified": 0,
        "supersedes": "R05/R06 H=18 (computed with the K~-line contraction artefact); R05's 18 representatives are not "
                      "in the corrected Nbar_g",
        "meaning": "invariant evanescent cohomology at O(kappa^3) on the complete local Y-space, formal H/T algebra, "
                   "parity-even; 'physical' only in this sense"})
    w("R07_H_basis", {"dim": len(rc["hatgen_H_basis_vectors"]), "basis_vectors": rc["hatgen_H_basis_vectors"]})
    w("R07_front1_reduced_hhc", {"R05rule": {k: h1[k] for k in h1 if k != "Rbar_sha256"}, "hatgen": f1,
        "Rbar_definition": "A + B2 + B3 + C + D2_EH + D3_EH (gauge-fixing part of S_2 and cbar channel E removed)",
        "outcome": "CLOSED" if (f1["solvable_with_ghost_eq"] and f1["verified_rows_zero"] and f1["particular_CT_all_evanescent"]
                                and h1["full_minus_reduced_equals_Dgf_plus_sigmaE"] and h1["tree"]["A_plus_B_EH_zero"]) else "OPEN"})
    idok = all(v.get("p22_CT_equals_minus_pole") for v in pz["identification"].values()) and \
        all(pz["identification"][g]["same_ratio"] for g in ("D01", "D03", "D11"))
    c = pz["completion"]
    w("R07_p22", {"identification": pz["identification"], "completion": c,
        "outcome": "CLOSED (three-point, one-loop pole level, reduced hat-aware system)" if
        (idok and c["solvable_with_ghost_eq"] and c["verified_rows_zero"] and c["completion_all_evanescent"]) else "OPEN",
        "ambiguity": f"the completion is fixed modulo Nbar_g; its class in Hbar (dim {hat['dim_H']}) is not determined by p22",
        "not_implied": "four-point rows, finite parts, two loops"})
    out1 = json.loads((CERT / "R07_front1_reduced_hhc.json").read_text())["outcome"]
    out2 = json.loads((CERT / "R07_p22.json").read_text())["outcome"]
    w("R07_P6_status", {"landmark": "XVI-R", "release": "R07", "P6_outcome": "Partial", "claim_status": "partial",
        "goals": {"1_reduced_hhc_loop_residual": out1, "2_p22_kappa2_consistency": out2.split(" ")[0],
                  "parity_odd_cohomology": "OPEN", "gram_reduced_scheme_count": "OPEN", "four_point_rows": "OPEN"},
        "O1_closure_unchanged": "CLOSED (scoped)", "bulk_symmetry_restored": False, "insertion": "frozen",
        "disclosed": ["R07-D1: hat-aware contraction extended to the external K~<->E graviton-probe line (and ghost line X<->Z); "
                      "reduced hhc linear block and sigma-bar images (F, KKc) rebuilt; Hbar 18 -> "
                      f"{hat['dim_H']} (T {hat['dim_T']}, blind {hat['dim_blind_g']}); R05/R06 H=18 superseded; frozen code "
                      "not edited",
                      "R07-D2: loop channels C and D_EH of the hhc row change under the extended rule (evanescent only; "
                      "4D closure unaffected)",
                      "R07-D3: reduced hhc loop residual Rbar defined as A+B+C+D_EH; full - reduced = D_gf + sigma_F E "
                      "(loop, purely evanescent)", "R07-D4: R07 code in r07/ (R01-R06 gates untouched)"]})
    M = {"NB": hat["dim_Nbar_g"], "TB": hat["dim_T"], "BG": hat["dim_blind_g"], "TCB": hat["dim_T_cap_blind"],
         "HB": hat["dim_H"], "RTRIV": hat["R05_18_trivial_count"], "BRo": hg["reduced_hhc_block"]["rank_old"],
         "BRn": hg["reduced_hhc_block"]["rank_new"], "BRj": hg["reduced_hhc_block"]["rank_joint"],
         "GFE": h1["Dgf_plus_sigmaE_n_terms"], "RBo": h1["Rbar_n_terms"], "RBn": f1["Rbar_n_terms"],
         "FKa": f1["particular_CT_n_terms"]["K"], "FCa": f1["particular_CT_n_terms"]["C"], "FHa": f1["particular_CT_n_terms"]["H"],
         "XIh": c["Xi_terms"]["hhc"], "XIk": c["Xi_terms"]["ccK"], "XIb": c["Xi_terms"]["ccb"],
         "PKa": c["completion_n_terms"]["K"], "PCa": c["completion_n_terms"]["C"], "PHa": c["completion_n_terms"]["H"],
         "NGF": hat["n_generators"]["F"], "NGG": hat["n_generators"]["G"], "NGK": hat["n_generators"]["KKc"],
         "PEC": hat["dim_pure_evanescent_combos"]}
    (R / "build/R07_numbers.tex").write_text("".join(f"\\newcommand{{\\{k}}}{{{v}}}\n" for k, v in M.items()))
    print("R07_CERTS_OK")


if __name__ == "__main__":
    main()
