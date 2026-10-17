#!/usr/bin/env python3
"""XVI-R R04 certificates + LaTeX fragments."""
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
    gd, gp, sh, fn = (j(x) for x in ("r04_gderive", "r04_gpat", "r04_s0hat", "r04_final"))
    w("R04_scstar_derivation", {"derivation": fn["derivation"], "kinetic_scan_nullspace": gd["coeff_nullspace"],
        "candidates_failing": {k: v.get("failing") for k, v in gd["candidates"].items()},
        "s0hat_tests": sh,
        "reading": "Derived sigma-bar on Y=c*.G from S-bar(Ktilde): X_K=(D^T Ktilde).G, X_C=rho*T(X_K), "
                   "X_S=-s0h.dY/dh (eps_S=-1). Raw kinetic cbar.Box.G is excluded (fails all 33). "
                   "s0 hat{h}=0 truncations of X_S make the count worse."})
    w("R04_G_outcome", {"G": fn["G"], "pattern_summary": {"n_ok": gp["n_ok"], "n_fail": gp["n_fail"],
        "fail_all_ccK_ccb": True, "all_failures_have_single_T_on_cstar": True},
        "outcome": "PARTIAL: derived image closes 20/33; 13/33 OPEN as T(c*) obstruction; no coefficient fix"})
    w("R04_classification", {"reduced": fn["reduced"], "transport_R03_18": fn["transport_18"],
        "H_bar_dim": fn["reduced"]["dim_H_bar"], "physical_certified": 0,
        "status": "CONDITIONAL on the T(c*) G-sector obstruction and on Y-space completeness",
        "note": "Failing G images excluded from the trivial span (their sigma-bar is not a BRST image); "
                "H_bar grows from R03's 18 to 20. R03's 18 all survive as nontrivial."})
    w("R04_H_bar_basis", {"dim": len(fn["H_bar_basis"]), "basis": fn["H_bar_basis"]})
    w("R04_P6_status", {"landmark": "XVI-R", "release": "R04", "P6_outcome": "Partial", "claim_status": "partial",
        "goals": {"1_scstar_derived": "ESTABLISHED: X_K, X_C=rho*T, X_S with eps_S=-1 from the reduced master",
                  "2_G_defect": "Partial: 20/33 closed; 13/33 OPEN (T(c*) obstruction on ghost rows)",
                  "3_classification": f"Partial: H_bar={fn['reduced']['dim_H_bar']} conditional classes "
                                      f"(R03's 18 survive + 2); 0 physical",
                  "4_Y_completeness": "OPEN (power-counting monomials only)",
                  "5_p22_kappa2_consistency": "OPEN"},
        "O1_closure_unchanged": "CLOSED (scoped)", "bulk_symmetry_restored": False, "insertion": "frozen",
        "disclosed": ["R04-D1: failing G images excluded from im(sigma-bar); H_bar 18->20 vs R03",
                      "R04-D2: R04 code in r04/ (R01-R03 gates untouched)"],
        "not_implied": "no direction certified physical; Missing never zero"})
    M = {"NP": fn["G"]["n_pass"], "NF": fn["G"]["n_fail"], "NGens": fn["G"]["n_images"],
         "HB": fn["reduced"]["dim_H_bar"], "TB": fn["reduced"]["dim_T_bar"], "BG": fn["reduced"]["dim_blind_g"],
         "TCB": fn["reduced"]["dim_T_cap_blind"], "NB": fn["reduced"]["dim_Nbar_g"],
         "TRS": fn["transport_18"]["nontrivial_in_H_bar"], "RS": fn["reduced"]["rank_Sigma_bar_passing"]}
    (R / "build/R04_numbers.tex").write_text("".join(f"\\newcommand{{\\{k}}}{{{v}}}\n" for k, v in M.items()))
    B = []
    for i, b in enumerate(fn["H_bar_basis"], 1):
        B.append(f"\\subsection*{{$\\bar{{\\mathcal O}}_{{{i}}}$: sectors {', '.join(b['sectors'])}; {b['n_terms']} terms}}")
        B.append(r"{\small\begin{align*}")
        B.append("\\\\\n".join(f"&{('+' if not c.startswith('-') else '')}{c}\\,[{s}]\\,{m}" for s, m, c in b["terms"]))
        B.append(r"\end{align*}}")
    (R / "build/R04_basis.tex").write_text("\n".join(B) + "\n")
    print("R04_CERTS_OK")


if __name__ == "__main__":
    main()
