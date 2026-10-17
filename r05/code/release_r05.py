#!/usr/bin/env python3
"""XVI-R R05 release gate. Layers: R03 ACCEPT (implies R02, R01, XVI-Q Q11)."""
import hashlib, json, sys
from pathlib import Path
R = Path(__file__).resolve().parents[2]
CERT = R / "certificates"
EXPECTED_SHA256 = "457e5f7cdf94267be31c7aa2094548471a265859f77acf6ef141702c4b834ab5"
SELF = "r05/code/release_r05.py"
NAMES = ("R05_obstruction_outcome", "R05_classification", "R05_H_basis", "R05_P6_status")
sys.path.insert(0, str(R / "r04/code"))
import release_r04 as G3


def sha(p):
    return hashlib.sha256((R / p).read_bytes()).hexdigest()


def tracked():
    out = set(G3.tracked()) | {"r04/code/release_r04.py", "r04/code/r04_expected.json", "r04/make_expected_r04.py"}
    out |= {p.relative_to(R).as_posix() for p in (R / "r05").rglob("*")
            if p.is_file() and p.suffix in (".py", ".tex") and "__pycache__" not in p.parts}
    return out - {SELF, "r05/code/r05_expected.json", "r05/make_expected_r05.py"} | {"r05/make_expected_r05.py"}


def checks():
    c = {}
    c["R04_layer_accept"] = json.loads((CERT / "R04_release.json").read_text())["decision"] == "ACCEPT"
    ep = R / "r05/code/r05_expected.json"
    c["expected_pinned"] = ep.is_file() and sha("r05/code/r05_expected.json") == EXPECTED_SHA256
    exp = json.loads(ep.read_text()) if ep.is_file() else {}
    certs = exp.get("certificates", {})
    c["certs_match_pinned"] = set(certs) == set(NAMES) and all(
        G3.G3.G2.strip(json.loads((CERT / f"{n}.json").read_text())) == v for n, v in certs.items())
    fs = exp.get("file_sha256", {})
    c["files_rehashed"] = bool(fs) and all((R / f).is_file() and sha(f) == h for f, h in fs.items())
    c["no_unpinned_file"] = tracked() <= set(fs)
    oo = json.loads((CERT / "R05_obstruction_outcome.json").read_text())
    c["closure_values"] = all(v["failing"] == 0 for v in oo["images"].values()) and oo["images"]["G"]["n"] == 33 and \
        all(oo["sanity_old_blocks_equal_Q11"].values()) and oo["outcome"].startswith("ARTEFACT")
    cl = json.loads((CERT / "R05_classification.json").read_text())
    c["classification_consistent"] = cl["H_dim"] == 18 and cl["physical_certified"] == 0 and \
        cl["reduced_hat"]["sigma_inside_Nbar"] is True and cl["transport_R04_20"]["nontrivial"] == 18 and \
        json.loads((CERT / "R05_H_basis.json").read_text())["dim"] == 18
    st = json.loads((CERT / "R05_P6_status.json").read_text())
    c["status_honest"] = st["P6_outcome"] == "Partial" and st["bulk_symmetry_restored"] is False and \
        st["insertion"] == "frozen" and st["goals"]["4_Y_completeness"] == "OPEN" and \
        st["goals"]["6_p22_kappa2_consistency"] == "OPEN"
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
    rec = {"landmark": "XVI-R", "release": "R05", "decision": "ACCEPT" if ok else "REJECT", "checks": c}
    (CERT / "R05_release.json").write_text(json.dumps(rec, indent=1) + "\n")
    print(json.dumps(rec, indent=1))
    if not ok:
        return 1
    (R / "RELEASE_STATUS.json").write_text(json.dumps({"landmark": "XVI-R", "release": "R05",
        "status": "released (honest partial)", "goals": st["goals"], "O1_closure": "CLOSED (scoped)",
        "bulk_symmetry_restored": False, "claim_status": "partial"}, indent=1) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
