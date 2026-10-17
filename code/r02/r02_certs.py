#!/usr/bin/env python3
"""XVI-R R02 certificates + generated LaTeX body (build/R02_body.tex)."""
import json
from pathlib import Path
R = Path(__file__).resolve().parents[2]
CERT = R / "certificates"


def strip(o):
    if isinstance(o, dict):
        return {k: strip(v) for k, v in o.items() if k != "seconds"}
    if isinstance(o, list):
        return [strip(v) for v in o]
    return o


def w(n, o):
    (CERT / f"{n}.json").write_text(json.dumps(o, indent=1, sort_keys=True, default=str) + "\n")


def main():
    s = strip(json.loads((R / "build/logs/r02_sci.json").read_text()))
    cl = s["classification_43"]
    w("R02_sigma_map", {k: s[k] for k in ("n_Y_generators", "rank_Sigma", "dim_kernel_constraints_on_Y",
                                          "dim_T_trivial_in_Ng", "dim_T_rows_only", "by_family", "leak")})
    w("R02_classification", {"dim_N_g": s["dim_N_g"], "dim_W_Q11": s["dim_W_Q11"], "unclassified_Q11": s["unclassified_Q11"],
      "dim_Wprime_cap_Ng": s["dim_Wprime_cap_Ng"], "dim_Ng_cap_spanSigma_plus_blind": s["dim_Ng_cap_spanSigma_plus_blind"],
      "classification_43": cl,
      "status": "CONDITIONAL: trivial = in the computed sigma-span satisfying all constraints; candidates are NOT proven "
                "nontrivial in cohomology because the Q11 composition gap (sigma(Y) not inside N) is open",
      "physical": "none certified"})
    w("R02_candidate_basis", {"dim": s["dim_H_candidate"], "robust": s["dim_H_candidate_robust"],
      "basis": s["basis"], "convention": "rref of N_g modulo W' in the Q08 coordinate order (canonical)"})
    w("R02_P6_status", {"landmark": "XVI-R", "release": "R02", "P6_outcome": "Partial", "claim_status": "partial",
      "goals": {"1_computed_trivial_span": "ESTABLISHED (computed sigma map on F/G/P monomial Y-space)",
                "2_classify_43": f"Partial: {cl['trivial_computed_sigma']} trivial (computed sigma), {cl['blind']} blind, "
                                 f"{cl['cohomology_candidate']} cohomology candidates; 0 physical",
                "3_invariant_basis": f"Partial: canonical {s['dim_H_candidate']}-dim candidate basis; invariance "
                                     "(nontriviality) conditional on closing the sigma-composition gap",
                "4_p22_kappa2_consistency": "OPEN"},
      "O1_closure_unchanged": "CLOSED (scoped)", "bulk_symmetry_restored": False, "insertion": "frozen",
      "disclosed_revisions": ["R02-R1: R01 gate's unpinned-file check excludes R02 files (pinned by R02 gate)"],
      "not_implied": "no direction certified physical; Missing never zero"})
    L = []
    L.append(r"\begin{center}\begin{tabular}{lrrr}\toprule family & generators & rank $\sigma$ & trivial in $N_g$\\\midrule")
    for k, v in s["by_family"].items():
        L.append(f"{k} & {v['n_generators']} & {v['rank_sigma']} & {v['dim_trivial_in_Ng']}\\\\")
    L.append(f"all & {s['n_Y_generators']} & {s['rank_Sigma']} & {s['dim_T_trivial_in_Ng']}\\\\\\bottomrule\\end{{tabular}}\\end{{center}}")
    (R / "build/R02_tab.tex").write_text("\n".join(L) + "\n")
    M = {"NY": s["n_Y_generators"], "RS": s["rank_Sigma"], "ZY": s["dim_kernel_constraints_on_Y"], "DT": s["dim_T_trivial_in_Ng"],
         "DTA": s["dim_T_rows_only"], "WQ": s["dim_W_Q11"], "WP": s["dim_Wprime_cap_Ng"], "WR": s["dim_Ng_cap_spanSigma_plus_blind"],
         "NT": cl["trivial_computed_sigma"], "NHC": cl["cohomology_candidate"], "NR": cl["of_which_robust"]}
    (R / "build/R02_numbers.tex").write_text("".join(f"\\newcommand{{\\{k}}}{{{v}}}\n" for k, v in M.items()))
    B = []
    for i, b in enumerate(s["basis"], 1):
        B.append(f"\\subsection*{{$\\mathcal O_{{{i}}}$: sectors {', '.join(b['sectors'])}; {b['n_terms']} terms}}")
        B.append(r"{\small\begin{align*}")
        rows = [f"&{('+' if not c.startswith('-') else '')}{c}\\,[{sec}]\\,{m}" for sec, m, c in b["terms"]]
        B.append("\\\\\n".join(rows))
        B.append(r"\end{align*}}")
    (R / "build/R02_basis.tex").write_text("\n".join(B) + "\n")
    print("R02_CERTS_OK")


if __name__ == "__main__":
    main()
