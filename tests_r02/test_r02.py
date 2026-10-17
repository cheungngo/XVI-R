import json
from pathlib import Path
R = Path(__file__).resolve().parents[1]


def j(n):
    return json.loads((R / "certificates" / f"{n}.json").read_text())


def test_counts_add_up():
    c = j("R02_classification")
    k = c["classification_43"]
    assert c["unclassified_Q11"] == 43 and k["trivial_computed_sigma"] + k["blind"] + k["cohomology_candidate"] == 43
    assert k["physical_proven"] == 0


def test_basis_dim():
    b = j("R02_candidate_basis")
    assert b["dim"] == len(b["basis"]) == j("R02_classification")["classification_43"]["cohomology_candidate"]


def test_no_leak():
    assert all(v == 0 for v in j("R02_sigma_map")["leak"].values())


def test_status():
    s = j("R02_P6_status")
    assert s["P6_outcome"] == "Partial" and s["bulk_symmetry_restored"] is False and s["insertion"] == "frozen"
