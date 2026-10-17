#!/usr/bin/env python3
"""XVI-R R01 certificates from build/logs/r01_sci.json and the XVI-Q Q11 layer rebuilt from the pinned zip."""
import hashlib, json, zipfile  # noqa: E401
from pathlib import Path
R = Path(__file__).resolve().parents[2]
QQ = R / "build" / "xviq" / "XVI-Q"
CERT = R / "certificates"
ZIP = R / "imports/XVI-Q-Q11-source-only.zip"


def strip(o):
    if isinstance(o, dict):
        return {k: strip(v) for k, v in o.items() if k != "seconds"}
    if isinstance(o, list):
        return [strip(v) for v in o]
    return o


def w(name, o):
    (CERT / f"{name}.json").write_text(json.dumps(o, indent=1, sort_keys=True, default=str) + "\n")


def main():
    CERT.mkdir(exist_ok=True)
    z = zipfile.ZipFile(ZIP)
    ref = json.loads(z.read("XVI-Q/code/p6q06/q07/q08/q09/q10/q11/q11_expected.json"))["certificates"]
    diff = []
    for n, o in ref.items():
        p = QQ / "certificates" / f"{n}.json"
        if not p.is_file():
            diff.append(n)
            continue
        a, b = strip(json.loads(p.read_text())), dict(o)
        if n in ("Q05_guard_scan", "Q09_guard_scan3"):
            a.pop("n_files_scanned", None)
            b.pop("n_files_scanned", None)
        if a != b:
            diff.append(n)
    q11r = json.loads((QQ / "certificates/Q11_release.json").read_text())
    w("R01_xviq_layer", {"xviq_source_zip_sha256": hashlib.sha256(ZIP.read_bytes()).hexdigest(),
                         "n_reference_certificates": len(ref), "mismatched_certificates": sorted(diff),
                         "Q11_release_decision_rebuilt": q11r["decision"],
                         "note": "XVI-Q Q11 rebuilt from the pinned source-only zip inside build/xviq; certificates compared "
                                 "to the values pinned in the zip's q11_expected.json (seconds stripped; guard-scan file counts excluded as in Q10-R1/Q11-R1)"})
    s = json.loads((R / "build/logs/r01_sci.json").read_text())
    g1 = s["goal1_defect_localisation"]
    w("R01_gap_localisation", {"results": strip(g1), "set_sizes": s["set_sizes"],
      "reading": "A one-sector correction delta with A(v+delta)=0 LOCALISES the hhc defect; it is not a derived sigma-image "
                 "and certifies nothing as trivial. Gap status: CHARACTERISED (localisation), NOT CLOSED."})
    g2 = s["goal2"]
    w("R01_unclassified_visibility", {**g2, "reading": "q(S) = number of unclassified N_g directions that cannot be "
      "represented with vanishing S-coordinates (modulo sigma-span + blind). Not a classification: no direction is "
      "certified trivial or physical."})
    w("R01_P6_status", {"landmark": "XVI-R", "release": "R01", "P6_outcome": "Partial", "claim_status": "partial",
      "goals": {"1_sigma_gap": "Partial: characterised by defect localisation; NOT CLOSED",
                "2_unclassified_43": "Partial: invariant sector/evanescent visibility computed; classification OPEN",
                "3_invariant_basis": "OPEN (not attempted beyond goal 2 visibility)",
                "4_p22_kappa2_consistency": "OPEN (not attempted)"},
      "O1_closure_unchanged": "CLOSED (scoped)", "bulk_symmetry_restored": False, "insertion": "frozen",
      "frozen_layer_revisions": [], "premises_imported": ["XVI-Q Q11 (pinned zip)"],
      "not_implied": "no direction certified physical or trivial; Missing never zero"})
    rows = []
    for fam in ("F", "G"):
        st = g1[fam]
        for k in ("sector_K", "sector_C", "sector_S", "sector_H", "CS", "hat_deg_1", "hat_deg_ge2", "all_coords"):
            if k in st:
                rows.append(f"{fam} & {k.replace('_', ' ')} & {st['n_fail']} & {st[k]['correctable']} & "
                            f"{st[k]['correctable_and_ghost_eq']} & {st[k]['defect_span_inside_image']} \\\\")
    vis = " ; ".join(f"{k.replace('_', ' ')}: {v}" for k, v in g2["visibility_q"].items())
    (R / "build/R01_numbers.tex").write_text("\\newcommand{\\GapRows}{" + "\n".join(rows) + "}\n"
        + "\\newcommand{\\Vis}{" + vis + "}\n" + f"\\newcommand{{\\DefRankF}}{{{g1['F'].get('rank_defect_span')}}}\n"
        + f"\\newcommand{{\\DefRankG}}{{{g1['G'].get('rank_defect_span')}}}\n")
    print("R01_CERTS_OK")


if __name__ == "__main__":
    main()
