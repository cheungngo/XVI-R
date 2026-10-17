#!/usr/bin/env python3
"""Freeze: pin R01 certificates and tracked files into code/r02/r02_expected.json and pin its sha in the gate."""
import hashlib, json, re, sys  # noqa: E401
from pathlib import Path
R = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(R / "code/r02"))
import release_r02 as G  # noqa: E402
certs = {n: G.strip(json.loads((R / f"certificates/{n}.json").read_text())) for n in
         ("R02_sigma_map", "R02_classification", "R02_candidate_basis", "R02_P6_status")}
fs = {f: G.sha(f) for f in sorted(G.tracked())}
p = R / "code/r02/r02_expected.json"
p.write_text(json.dumps({"certificates": certs, "file_sha256": fs}, indent=1, sort_keys=True) + "\n")
h = hashlib.sha256(p.read_bytes()).hexdigest()
g = R / G.SELF
g.write_text(re.sub(r'EXPECTED_SHA256 = "[^"]*"', f'EXPECTED_SHA256 = "{h}"', g.read_text()))
print(h)
