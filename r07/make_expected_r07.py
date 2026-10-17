#!/usr/bin/env python3
import hashlib, json, re, sys
from pathlib import Path
R = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(R / "r07/code"))
import release_r07 as G
certs = {n: G.G3.G3.G3.G3.G2.strip(json.loads((R / f"certificates/{n}.json").read_text())) for n in G.NAMES}
fs = {f: G.sha(f) for f in sorted(G.tracked() - {"r07/code/r07_expected.json"})}
p = R / "r07/code/r07_expected.json"
p.write_text(json.dumps({"certificates": certs, "file_sha256": fs}, indent=1, sort_keys=True) + "\n")
h = hashlib.sha256(p.read_bytes()).hexdigest()
g = R / G.SELF
g.write_text(re.sub(r'EXPECTED_SHA256 = "[^"]*"', f'EXPECTED_SHA256 = "{h}"', g.read_text()))
print(h)
