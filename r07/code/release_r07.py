#!/usr/bin/env python3
"""XVI-R R07 release gate. Layers: R06 ACCEPT (implies R05..R01, XVI-Q Q11)."""
import hashlib, json, sys
from pathlib import Path
R = Path(__file__).resolve().parents[2]
CERT = R / "certificates"
EXPECTED_SHA256 = "b0dfd8850c7d39b5678890da2932c2ac3aedd29324262891bc2b16fc994fc8a3"
SELF = "r07/code/release_r07.py"
NAMES = ("R07_contraction_Kline", "R07_classification", "R07_H_basis", "R07_front1_reduced_hhc", "R07_p22", "R07_P6_status")
sys.path.insert(0, str(R / "r06/code"))
import release_r06 as G3


def sha(p):
    return hashlib.sha256((R / p).read_bytes()).hexdigest()


def tracked():
    out = set(G3.tracked()) | {"r06/code/release_r06.py", "r06/code/r06_expected.json", "r06/make_expected_r06.py"}
    out |= {p.relative_to(R).as_posix() for p in (R / "r07").rglob("*")
            if p.is_file() and p.suffix in (".py", ".tex") and "__pycache__" not in p.parts}
    return out - {SELF, "r07/code/r07_expected.json", "r07/make_expected_r07.py"} | {"r07/make_expected_r07.py"}


def checks():
    c = {}
    c["R06_layer_accept"] = json.loads((CERT / "R06_release.json").read_text())["decision"] == "ACCEPT"
    ep = R / "r07/code/r07_expected.json"
    c["expected_pinned"] = ep.is_file() and sha("r07/code/r07_expected.json") == EXPECTED_SHA256
    exp = json.loads(ep.read_text()) if ep.is_file() else {}
    certs = exp.get("certificates", {})
    c["certs_match_pinned"] = set(certs) == set(NAMES) and all(
        G3.G3.G3.G3.G2.strip(json.loads((CERT / f"{n}.json").read_text())) == v for n, v in certs.items())
    fs = exp.get("file_sha256", {})
    c["files_rehashed"] = bool(fs) and all((R / f).is_file() and sha(f) == h for f, h in fs.items())
    c["no_unpinned_file"] = tracked() <= set(fs)
    kl = json.loads((CERT / "R07_contraction_Kline.json").read_text())
    c["Kline_finding_recorded"] = kl["reduced_hhc_block"]["rowspace_changed"] is True and kl["loop_chan_C_change_4D"] == "0" \
        and kl["loop_chan_D_EH_change_4D"] == "0"
    cl = json.loads((CERT / "R07_classification.json").read_text())
    c["classification_consistent"] = cl["R05rule_sanity"]["dim_H"] == 18 and cl["hatgen"]["sigma_inside_Nbar"] is True and \
        cl["hatgen"]["images_inside_Nbar"] is True and cl["H_dim"] == cl["H_basis_dim"] == \
        json.loads((CERT / "R07_H_basis.json").read_text())["dim"] and cl["physical_certified"] == 0
    f1 = json.loads((CERT / "R07_front1_reduced_hhc.json").read_text())
    c["front1_consistent"] = f1["outcome"] in ("CLOSED", "OPEN") and (f1["outcome"] == "OPEN" or (
        f1["hatgen"]["solvable_with_ghost_eq"] and f1["hatgen"]["verified_rows_zero"]))
    pz = json.loads((CERT / "R07_p22.json").read_text())
    c["p22_consistent"] = pz["outcome"].startswith("OPEN") or (pz["completion"]["solvable_with_ghost_eq"] and
        pz["completion"]["verified_rows_zero"] and pz["completion"]["completion_all_evanescent"])
    st = json.loads((CERT / "R07_P6_status.json").read_text())
    c["status_honest"] = st["P6_outcome"] == "Partial" and st["bulk_symmetry_restored"] is False and \
        st["insertion"] == "frozen" and st["goals"]["parity_odd_cohomology"] == "OPEN"
    reg = json.loads((R / "registration/local_registration.json").read_text())["hard_stops"]
    c["hard_stops"] = reg["bulk_symmetry_restored"] is False and reg["no_workspace_sync"] is True
    return c, st


def main():
    (R / "RELEASE_STATUS.json").unlink(missing_ok=True)
    try:
        c, st = checks()
    except Exception as e:
        c, st = {"exception": False, "error": repr(e)}, {}
    ok = all(v is True for k, v in c.items() if k != "error")
    rec = {"landmark": "XVI-R", "release": "R07", "decision": "ACCEPT" if ok else "REJECT", "checks": c}
    (CERT / "R07_release.json").write_text(json.dumps(rec, indent=1) + "\n")
    print(json.dumps(rec, indent=1))
    if not ok:
        return 1
    (R / "RELEASE_STATUS.json").write_text(json.dumps({"landmark": "XVI-R", "release": "R07",
        "status": "released (honest partial)", "H": f"{json.loads((CERT / 'R07_classification.json').read_text())['H_dim']} = invariant evanescent cohomology at O(kappa^3) (corrected, R07-D1)", "goals": st["goals"], "O1_closure": "CLOSED (scoped)",
        "bulk_symmetry_restored": False, "claim_status": "partial"}, indent=1) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
