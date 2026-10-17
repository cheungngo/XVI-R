#!/usr/bin/env python3
"""XVI-R R08 release gate. Layers: R07 ACCEPT (implies R06..R01, XVI-Q Q11)."""
import hashlib, json, sys
from pathlib import Path
R = Path(__file__).resolve().parents[2]
CERT = R / "certificates"
EXPECTED_SHA256 = "da91f5d0ceb4374c1bf6f60384aa0f5b21219bf0fde9eb9892bccdf59e420889"
SELF = "r08/code/release_r08.py"
NAMES = ("R08_gram", "R08_parity", "R08_four_point", "R08_P6_status")
sys.path.insert(0, str(R / "r07/code"))
import release_r07 as G3


def sha(p):
    return hashlib.sha256((R / p).read_bytes()).hexdigest()


def tracked():
    out = set(G3.tracked()) | {"r07/code/release_r07.py", "r07/code/r07_expected.json", "r07/make_expected_r07.py"}
    out |= {p.relative_to(R).as_posix() for p in (R / "r08").rglob("*")
            if p.is_file() and p.suffix in (".py", ".tex") and "__pycache__" not in p.parts}
    return out - {SELF, "r08/code/r08_expected.json", "r08/make_expected_r08.py"} | {"r08/make_expected_r08.py"}


def checks():
    c = {}
    c["R07_layer_accept"] = json.loads((CERT / "R07_release.json").read_text())["decision"] == "ACCEPT"
    ep = R / "r08/code/r08_expected.json"
    c["expected_pinned"] = ep.is_file() and sha("r08/code/r08_expected.json") == EXPECTED_SHA256
    exp = json.loads(ep.read_text()) if ep.is_file() else {}
    certs = exp.get("certificates", {})
    c["certs_match_pinned"] = set(certs) == set(NAMES) and all(
        G3.G3.G3.G3.G3.G2.strip(json.loads((CERT / f"{n}.json").read_text())) == v for n, v in certs.items())
    fs = exp.get("file_sha256", {})
    c["files_rehashed"] = bool(fs) and all((R / f).is_file() and sha(f) == h for f, h in fs.items())
    c["no_unpinned_file"] = tracked() <= set(fs)
    g = json.loads((CERT / "R08_gram.json").read_text())
    c["gram_consistent"] = g["R07_reproduced"]["dim_H"] == 6 and g["G_sigma_image_zero"] is True and \
        g["H_gram_reduced"] == g["gram_scheme"]["dim_H"] and g["outcome"] in ("CLOSED", "OPEN")
    p = json.loads((CERT / "R08_parity.json").read_text())
    c["parity_honest"] = p["odd_cohomology_dimension"].startswith("Missing") and p["outcome"] in ("PARTIAL", "OPEN")
    f = json.loads((CERT / "R08_four_point.json").read_text())
    c["four_point_honest"] = f["outcome"].startswith("OPEN") and len(f["missing"]) > 0
    st = json.loads((CERT / "R08_P6_status.json").read_text())
    c["status_honest"] = st["P6_outcome"] == "Partial" and st["bulk_symmetry_restored"] is False and \
        st["insertion"] == "frozen" and st["goals"]["1_four_point_rows"] == "OPEN"
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
    rec = {"landmark": "XVI-R", "release": "R08", "decision": "ACCEPT" if ok else "REJECT", "checks": c}
    (CERT / "R08_release.json").write_text(json.dumps(rec, indent=1) + "\n")
    print(json.dumps(rec, indent=1))
    if not ok:
        return 1
    (R / "RELEASE_STATUS.json").write_text(json.dumps({"landmark": "XVI-R", "release": "R08",
        "status": "released (honest partial)", "H": "6 (formal and Gram-reduced, parity-even); parity-odd Missing", "goals": st["goals"], "O1_closure": "CLOSED (scoped)",
        "bulk_symmetry_restored": False, "claim_status": "partial"}, indent=1) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
