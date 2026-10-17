#!/usr/bin/env python3
"""XVI-R R01 release gate (fails closed; on reject no RELEASE_STATUS.json is left).
Layer gate: XVI-Q Q11 rebuilt from the pinned zip, its own Q11 gate ACCEPTs and every certificate equals the zip's.
R01 gate: expected file pinned by sha; every R01 certificate equals its pinned copy; every code/tools/tests/
registration/intake/imports/manuscript-source file pinned and re-hashed (an unpinned file rejects); status honest."""
import hashlib, json  # noqa: E401
from pathlib import Path
R = Path(__file__).resolve().parents[2]
CERT = R / "certificates"
EXPECTED_SHA256 = "0d2857876177187021595fb41a9cff917741893b5a485e3cfbb57c36fedb13be"
SELF = "code/r01/release_r01.py"
ZIPSHA = "1888f3277b793ee8cbf3954b3819bec03e41d8c4dc59757744dfda8b147a35ad"


def sha(p):
    return hashlib.sha256((R / p).read_bytes()).hexdigest()


def strip(o):
    if isinstance(o, dict):
        return {k: strip(v) for k, v in o.items() if k != "seconds"}
    if isinstance(o, list):
        return [strip(v) for v in o]
    return o


def tracked():
    out = set()
    for d, pat in (("code", "**/*.py"), ("code", "**/*.json"), ("tools", "*.py"), ("tests", "*.py"), ("registration", "*"),
                   ("intake", "*"), ("imports", "*"), ("manuscript", "*.tex")):
        out |= {p.relative_to(R).as_posix() for p in (R / d).glob(pat) if p.is_file() and "__pycache__" not in p.parts}
    return out - {SELF, "code/r01/r01_expected.json"}


def checks():
    c = {}
    ep = R / "code/r01/r01_expected.json"
    c["expected_pinned"] = ep.is_file() and sha("code/r01/r01_expected.json") == EXPECTED_SHA256
    exp = json.loads(ep.read_text()) if ep.is_file() else {}
    certs = exp.get("certificates", {})
    c["certs_match_pinned"] = len(certs) == 4 and all(
        (CERT / f"{n}.json").is_file() and strip(json.loads((CERT / f"{n}.json").read_text())) == v for n, v in certs.items())
    fs = exp.get("file_sha256", {})
    c["files_rehashed"] = bool(fs) and all((R / f).is_file() and sha(f) == h for f, h in fs.items())
    # R02-R1 (disclosed): files introduced by R02 (path contains r02/R02) are pinned by the R02 gate, not here
    c["no_unpinned_file"] = {f for f in tracked() if "r02" not in f.lower()} <= set(fs)
    lay = json.loads((CERT / "R01_xviq_layer.json").read_text())
    c["xviq_layer"] = lay["xviq_source_zip_sha256"] == ZIPSHA and lay["mismatched_certificates"] == [] and \
        lay["Q11_release_decision_rebuilt"] == "ACCEPT" and lay["n_reference_certificates"] >= 44
    st = json.loads((CERT / "R01_P6_status.json").read_text())
    c["status_honest"] = st["P6_outcome"] == "Partial" and st["bulk_symmetry_restored"] is False and \
        st["insertion"] == "frozen" and st["O1_closure_unchanged"] == "CLOSED (scoped)" and \
        "NOT CLOSED" in st["goals"]["1_sigma_gap"] and st["goals"]["3_invariant_basis"].startswith("OPEN")
    reg = json.loads((R / "registration/local_registration.json").read_text())["hard_stops"]
    c["hard_stops"] = reg["bulk_symmetry_restored"] is False and reg["no_workspace_sync"] is True and \
        reg["Route_B_discharged"] is False and reg["XV_B_opened"] is False and reg["XIV_J_opened"] is False
    return c, st


def main():
    for p in (R / "RELEASE_STATUS.json",):
        p.unlink(missing_ok=True)
    try:
        c, st = checks()
    except Exception as e:  # noqa: BLE001
        c, st = {"exception": False, "error": repr(e)}, {}
    ok = all(v is True for k, v in c.items() if k != "error")
    rec = {"landmark": "XVI-R", "release": "R01", "decision": "ACCEPT" if ok else "REJECT", "checks": c}
    (CERT / "R01_release.json").write_text(json.dumps(rec, indent=1) + "\n")
    print(json.dumps(rec, indent=1))
    if not ok:
        return 1
    (R / "RELEASE_STATUS.json").write_text(json.dumps({"landmark": "XVI-R", "release": "R01", "status": "released (honest partial)",
        "goals": st["goals"], "O1_closure": "CLOSED (scoped)", "bulk_symmetry_restored": False, "claim_status": "partial"}, indent=1) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
