import json
from pathlib import Path
R = Path(__file__).resolve().parents[2]


def j(n):
    return json.loads((R / "certificates" / f"{n}.json").read_text())


def test_all_families_closed():
    o = j("R05_obstruction_outcome")
    assert all(v["failing"] == 0 for v in o["images"].values())
    assert o["G_eps_plus1_failing"] > 0


def test_sanity_frozen_blocks():
    assert all(j("R05_obstruction_outcome")["sanity_old_blocks_equal_Q11"].values())


def test_classification():
    c = j("R05_classification")
    assert c["H_dim"] == 18 and c["reduced_hat"]["sigma_inside_Nbar"] is True and c["physical_certified"] == 0
    t = c["transport_R04_20"]
    assert t["trivial_or_blind"] + t["nontrivial"] + t["not_cocycles"] == 20


def test_status():
    s = j("R05_P6_status")
    assert s["P6_outcome"] == "Partial" and s["bulk_symmetry_restored"] is False
