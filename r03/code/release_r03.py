#!/usr/bin/env python3
"""XVI-R R03 release gate (fails closed).  Layers: R02 gate ACCEPT (which requires R01 ACCEPT incl. the XVI-Q Q11 layer).
R03: expected file pinned; certificates equal pins; every tracked file (R01/R02 sets + r03/**) pinned and re-hashed;
values consistent; status honest (gap: hhc closed only in reduced variables, G-sector open; 0 physical)."""
import hashlib, json, sys  # noqa: E401
from pathlib import Path
R = Path(__file__).resolve().parents[2]
CERT = R / "certificates"
EXPECTED_SHA256 = "c01a86f1f23561d354af8df72507fd2b0bef924f311e76b232d217182a7386f2"
SELF = "r03/code/release_r03.py"
NAMES = ("R03_gap_diagnosis", "R03_reduced", "R03_classification", "R03_H_bar_basis", "R03_P6_status")
sys.path.insert(0, str(R / "code/r02"))
import release_r02 as G2  # noqa: E402


def sha(p):
    return hashlib.sha256((R / p).read_bytes()).hexdigest()


def tracked():
    out = set(G2.tracked()) | {"code/r02/release_r02.py", "code/r01/r01_expected.json", "code/r02/r02_expected.json"}
    out |= {p.relative_to(R).as_posix() for p in (R / "r03").rglob("*") if p.is_file() and p.suffix in (".py", ".tex")
            and "__pycache__" not in p.parts}
    return out - {SELF, "r03/code/r03_expected.json", "r03/make_expected_r03.py"} | {"r03/make_expected_r03.py"}


def checks():
    c = {}
    c["R02_layer_accept"] = json.loads((CERT / "R02_release.json").read_text())["decision"] == "ACCEPT"
    ep = R / "r03/code/r03_expected.json"
    c["expected_pinned"] = ep.is_file() and sha("r03/code/r03_expected.json") == EXPECTED_SHA256
    exp = json.loads(ep.read_text()) if ep.is_file() else {}
    certs = exp.get("certificates", {})
    c["certs_match_pinned"] = set(certs) == set(NAMES) and all(
        G2.strip(json.loads((CERT / f"{n}.json").read_text())) == v for n, v in certs.items())
    fs = exp.get("file_sha256", {})
    c["files_rehashed"] = bool(fs) and all((R / f).is_file() and sha(f) == h for f, h in fs.items())
    c["no_unpinned_file"] = tracked() <= set(fs)
    rd = json.loads((CERT / "R03_reduced.json").read_text())
    c["gap_values"] = rd["F_images_failing_reduced"] == 0 and rd["KKc"]["leak"] == 0 and \
        set(rd["KKc"]["blocks"]) == {"pass"} and "OPEN" in rd["gap_outcome"] and rd["dim_N_g"] == 104
    cl = json.loads((CERT / "R03_classification.json").read_text())
    t = cl["transport_R02_22"]
    c["classification_consistent"] = t["dim"] == 22 and t["fail_reduced_constraints"] + t["trivial_or_blind_in_reduced"] + \
        t["nontrivial_classes_in_H_bar"] == 22 and cl["physical_certified"] == 0 and cl["status"].startswith("CONDITIONAL") and \
        json.loads((CERT / "R03_H_bar_basis.json").read_text())["dim"] == cl["H_bar_dim_nonblind"]
    st = json.loads((CERT / "R03_P6_status.json").read_text())
    c["status_honest"] = st["P6_outcome"] == "Partial" and st["bulk_symmetry_restored"] is False and \
        st["insertion"] == "frozen" and st["O1_closure_unchanged"] == "CLOSED (scoped)" and \
        "OPEN" in st["goals"]["2_gap"] and st["goals"]["4_p22_kappa2_consistency"] == "OPEN"
    reg = json.loads((R / "registration/local_registration.json").read_text())["hard_stops"]
    c["hard_stops"] = reg["bulk_symmetry_restored"] is False and reg["no_workspace_sync"] is True
    return c, st


def main():
    (R / "RELEASE_STATUS.json").unlink(missing_ok=True)
    try:
        c, st = checks()
    except Exception as e:  # noqa: BLE001
        c, st = {"exception": False, "error": repr(e)}, {}
    ok = all(v is True for k, v in c.items() if k != "error")
    rec = {"landmark": "XVI-R", "release": "R03", "decision": "ACCEPT" if ok else "REJECT", "checks": c}
    (CERT / "R03_release.json").write_text(json.dumps(rec, indent=1) + "\n")
    print(json.dumps(rec, indent=1))
    if not ok:
        return 1
    (R / "RELEASE_STATUS.json").write_text(json.dumps({"landmark": "XVI-R", "release": "R03", "status": "released (honest partial)",
        "goals": st["goals"], "O1_closure": "CLOSED (scoped)", "bulk_symmetry_restored": False, "claim_status": "partial"},
        indent=1) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
