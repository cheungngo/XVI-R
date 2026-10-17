#!/usr/bin/env python3
"""XVI-R R04 release gate. Layers: R03 ACCEPT (implies R02, R01, XVI-Q Q11)."""
import hashlib, json, sys
from pathlib import Path
R = Path(__file__).resolve().parents[2]
CERT = R / "certificates"
EXPECTED_SHA256 = "5e2625eb12de66dc478bedb9ac8c2db6b46817ce937ae156f42d7c3211040139"
SELF = "r04/code/release_r04.py"
NAMES = ("R04_scstar_derivation", "R04_G_outcome", "R04_classification", "R04_H_bar_basis", "R04_P6_status")
sys.path.insert(0, str(R / "r03/code"))
import release_r03 as G3


def sha(p):
    return hashlib.sha256((R / p).read_bytes()).hexdigest()


def tracked():
    out = set(G3.tracked()) | {"r03/code/release_r03.py", "r03/code/r03_expected.json", "r03/make_expected_r03.py"}
    out |= {p.relative_to(R).as_posix() for p in (R / "r04").rglob("*")
            if p.is_file() and p.suffix in (".py", ".tex") and "__pycache__" not in p.parts}
    return out - {SELF, "r04/code/r04_expected.json", "r04/make_expected_r04.py"} | {"r04/make_expected_r04.py"}


def checks():
    c = {}
    c["R03_layer_accept"] = json.loads((CERT / "R03_release.json").read_text())["decision"] == "ACCEPT"
    ep = R / "r04/code/r04_expected.json"
    c["expected_pinned"] = ep.is_file() and sha("r04/code/r04_expected.json") == EXPECTED_SHA256
    exp = json.loads(ep.read_text()) if ep.is_file() else {}
    certs = exp.get("certificates", {})
    c["certs_match_pinned"] = set(certs) == set(NAMES) and all(
        G3.G2.strip(json.loads((CERT / f"{n}.json").read_text())) == v for n, v in certs.items())
    fs = exp.get("file_sha256", {})
    c["files_rehashed"] = bool(fs) and all((R / f).is_file() and sha(f) == h for f, h in fs.items())
    c["no_unpinned_file"] = tracked() <= set(fs)
    go = json.loads((CERT / "R04_G_outcome.json").read_text())
    c["G_values"] = go["G"]["n_pass"] == 20 and go["G"]["n_fail"] == 13 and "OPEN" in go["outcome"]
    cl = json.loads((CERT / "R04_classification.json").read_text())
    c["classification_consistent"] = cl["H_bar_dim"] == 20 and cl["physical_certified"] == 0 and \
        cl["transport_R03_18"]["nontrivial_in_H_bar"] == 18 and cl["status"].startswith("CONDITIONAL") and \
        json.loads((CERT / "R04_H_bar_basis.json").read_text())["dim"] == 20
    st = json.loads((CERT / "R04_P6_status.json").read_text())
    c["status_honest"] = st["P6_outcome"] == "Partial" and st["bulk_symmetry_restored"] is False and \
        st["insertion"] == "frozen" and "OPEN" in st["goals"]["2_G_defect"] and \
        st["goals"]["5_p22_kappa2_consistency"] == "OPEN"
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
    rec = {"landmark": "XVI-R", "release": "R04", "decision": "ACCEPT" if ok else "REJECT", "checks": c}
    (CERT / "R04_release.json").write_text(json.dumps(rec, indent=1) + "\n")
    print(json.dumps(rec, indent=1))
    if not ok:
        return 1
    (R / "RELEASE_STATUS.json").write_text(json.dumps({"landmark": "XVI-R", "release": "R04",
        "status": "released (honest partial)", "goals": st["goals"], "O1_closure": "CLOSED (scoped)",
        "bulk_symmetry_restored": False, "claim_status": "partial"}, indent=1) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
