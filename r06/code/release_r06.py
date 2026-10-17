#!/usr/bin/env python3
"""XVI-R R06 release gate. Layers: R05 ACCEPT (implies R04..R01, XVI-Q Q11)."""
import hashlib, json, sys
from pathlib import Path
R = Path(__file__).resolve().parents[2]
CERT = R / "certificates"
EXPECTED_SHA256 = "da0a940704976e6940c8e8831312881dbd7e253d77d893b38cd109c535d4ed00"
SELF = "r06/code/release_r06.py"
NAMES = ("R06_loop_remainders", "R06_completeness", "R06_classification", "R06_P6_status")
sys.path.insert(0, str(R / "r05/code"))
import release_r05 as G3


def sha(p):
    return hashlib.sha256((R / p).read_bytes()).hexdigest()


def tracked():
    out = set(G3.tracked()) | {"r05/code/release_r05.py", "r05/code/r05_expected.json", "r05/make_expected_r05.py"}
    out |= {p.relative_to(R).as_posix() for p in (R / "r06").rglob("*")
            if p.is_file() and p.suffix in (".py", ".tex") and "__pycache__" not in p.parts}
    return out - {SELF, "r06/code/r06_expected.json", "r06/make_expected_r06.py"} | {"r06/make_expected_r06.py"}


def checks():
    c = {}
    c["R05_layer_accept"] = json.loads((CERT / "R05_release.json").read_text())["decision"] == "ACCEPT"
    ep = R / "r06/code/r06_expected.json"
    c["expected_pinned"] = ep.is_file() and sha("r06/code/r06_expected.json") == EXPECTED_SHA256
    exp = json.loads(ep.read_text()) if ep.is_file() else {}
    certs = exp.get("certificates", {})
    c["certs_match_pinned"] = set(certs) == set(NAMES) and all(
        G3.G3.G3.G2.strip(json.loads((CERT / f"{n}.json").read_text())) == v for n, v in certs.items())
    fs = exp.get("file_sha256", {})
    c["files_rehashed"] = bool(fs) and all((R / f).is_file() and sha(f) == h for f, h in fs.items())
    c["no_unpinned_file"] = tracked() <= set(fs)
    lr = json.loads((CERT / "R06_loop_remainders.json").read_text())
    c["A_propagated"] = lr["loop_remainders_unchanged"] is False and lr["difference_purely_evanescent"] is True and \
        lr["refit"]["Q_full_hhc_hat"]["solvable_with_ghost_eq"] is True and lr["refit"]["Q_full_hhc_hat"]["verified_rows_zero"] is True
    cp = json.loads((CERT / "R06_completeness.json").read_text())
    c["B_complete"] = cp["leak_outside_ev_plus_4d"] == 0 and cp["transport_R05_18"]["stable"] is True and \
        cp["enlarged"]["new_trivial_directions"] == 0 and all(v["partition_ok"] and v["frozen_equals_independent_ev"] for v in cp["counts"].values())
    cl = json.loads((CERT / "R06_classification.json").read_text())
    c["classification_consistent"] = cl["H_dim"] == 18 and cl["physical_certified"] == 0 and cl["invariant_evanescent_cohomology"] is True
    st = json.loads((CERT / "R06_P6_status.json").read_text())
    c["status_honest"] = st["P6_outcome"] == "Partial" and st["bulk_symmetry_restored"] is False and \
        st["insertion"] == "frozen" and st["goals"]["p22_kappa2_consistency"] == "OPEN" and \
        st["goals"]["reduced_formulation_loop_residual"].startswith("OPEN")
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
    rec = {"landmark": "XVI-R", "release": "R06", "decision": "ACCEPT" if ok else "REJECT", "checks": c}
    (CERT / "R06_release.json").write_text(json.dumps(rec, indent=1) + "\n")
    print(json.dumps(rec, indent=1))
    if not ok:
        return 1
    (R / "RELEASE_STATUS.json").write_text(json.dumps({"landmark": "XVI-R", "release": "R06",
        "status": "released (honest partial)", "H": "18 = invariant evanescent cohomology at O(kappa^3)", "goals": st["goals"], "O1_closure": "CLOSED (scoped)",
        "bulk_symmetry_restored": False, "claim_status": "partial"}, indent=1) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
