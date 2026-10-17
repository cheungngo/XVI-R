import json
from pathlib import Path
R = Path(__file__).resolve().parents[2]


def j(n):
    return json.loads((R / "certificates" / f"{n}.json").read_text())


def test_derived_closes_20():
    g = j("R04_G_outcome")["G"]
    assert g["n_pass"] == 20 and g["n_fail"] == 13


def test_kinetic_ruled_out():
    d = j("R04_scstar_derivation")
    assert d["kinetic_scan_nullspace"] == []
    assert d["candidates_failing"]["Q08_handmade"] == 33


def test_classification():
    c = j("R04_classification")
    assert c["H_bar_dim"] == 20 and c["physical_certified"] == 0
    assert c["transport_R03_18"]["nontrivial_in_H_bar"] == 18


def test_status():
    s = j("R04_P6_status")
    assert s["P6_outcome"] == "Partial" and s["bulk_symmetry_restored"] is False
