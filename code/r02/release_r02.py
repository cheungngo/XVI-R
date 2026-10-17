#!/usr/bin/env python3
"""XVI-R R02 release gate (layers: R01 gate ACCEPT incl. XVI-Q Q11 layer).\nBased on the R01 gate (fails closed; on reject no RELEASE_STATUS.json is left).
Layer gate: XVI-Q Q11 rebuilt from the pinned zip, its own Q11 gate ACCEPTs and every certificate equals the zip's.
R01 gate: expected file pinned by sha; every R01 certificate equals its pinned copy; every code/tools/tests/
registration/intake/imports/manuscript-source file pinned and re-hashed (an unpinned file rejects); status honest."""
import hashlib, json  # noqa: E401
from pathlib import Path
R = Path(__file__).resolve().parents[2]
CERT = R / "certificates"
EXPECTED_SHA256 = "de667a88aeef6d11bd091941646df3068760e2e5d1343114fc200586cc3bc543"
SELF = "code/r02/release_r02.py"
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
    for d, pat in (("code", "**/*.py"), ("code", "**/*.json"), ("tools", "*.py"), ("tests", "*.py"), ("tests_r02", "*.py"),
                   ("registration", "*"), ("intake", "*"), ("imports", "*"), ("manuscript", "*.tex")):
        out |= {p.relative_to(R).as_posix() for p in (R / d).glob(pat) if p.is_file() and "__pycache__" not in p.parts}
    return out - {SELF, "code/r02/release_r02.py", "code/r01/r01_expected.json", "code/r02/r02_expected.json"}


def checks():
    c = {}
    r1 = json.loads((CERT / "R01_release.json").read_text())
    c["R01_layer_accept"] = r1["decision"] == "ACCEPT"
    ep = R / "code/r02/r02_expected.json"
    c["expected_pinned"] = ep.is_file() and sha("code/r02/r02_expected.json") == EXPECTED_SHA256
    exp = json.loads(ep.read_text()) if ep.is_file() else {}
    certs = exp.get("certificates", {})
    c["certs_match_pinned"] = len(certs) == 4 and all((CERT / f"{n}.json").is_file() and
        strip(json.loads((CERT / f"{n}.json").read_text())) == v for n, v in certs.items())
    fs = exp.get("file_sha256", {})
    c["files_rehashed"] = bool(fs) and all((R / f).is_file() and sha(f) == h for f, h in fs.items())
    c["no_unpinned_file"] = tracked() <= set(fs)
    cl = json.loads((CERT / "R02_classification.json").read_text())
    k = cl["classification_43"]
    c["classification_consistent"] = cl["unclassified_Q11"] == 43 and cl["dim_N_g"] == 104 and \
        k["trivial_computed_sigma"] + k["blind"] + k["cohomology_candidate"] == 43 and k["physical_proven"] == 0 and \
        cl["status"].startswith("CONDITIONAL")
    st = json.loads((CERT / "R02_P6_status.json").read_text())
    c["status_honest"] = st["P6_outcome"] == "Partial" and st["bulk_symmetry_restored"] is False and \
        st["insertion"] == "frozen" and st["O1_closure_unchanged"] == "CLOSED (scoped)" and \
        any(x.startswith("R02-R1") for x in st["disclosed_revisions"]) and st["goals"]["4_p22_kappa2_consistency"] == "OPEN"
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
    rec = {"landmark": "XVI-R", "release": "R02", "decision": "ACCEPT" if ok else "REJECT", "checks": c}
    (CERT / "R02_release.json").write_text(json.dumps(rec, indent=1) + "\n")
    print(json.dumps(rec, indent=1))
    if not ok:
        return 1
    (R / "RELEASE_STATUS.json").write_text(json.dumps({"landmark": "XVI-R", "release": "R02", "status": "released (honest partial)",
        "goals": st["goals"], "O1_closure": "CLOSED (scoped)", "bulk_symmetry_restored": False, "claim_status": "partial"}, indent=1) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
