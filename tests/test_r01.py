import json
from pathlib import Path
R = Path(__file__).resolve().parents[1]


def j(n):
    return json.loads((R / "certificates" / f"{n}.json").read_text())


def test_layer():
    a = j("R01_xviq_layer")
    assert a["mismatched_certificates"] == [] and a["Q11_release_decision_rebuilt"] == "ACCEPT"


def test_unclassified_consistent_with_q11():
    g = j("R01_unclassified_visibility")
    assert g["dim_N_g"] == 104 and g["unclassified"] == 43
    assert all(0 <= v <= 43 for v in g["visibility_q"].values())


def test_gap_counts_match_q11():
    g = j("R01_gap_localisation")["results"]
    assert g["F"]["n_fail"] == 123


def test_status():
    s = j("R01_P6_status")
    assert s["P6_outcome"] == "Partial" and s["bulk_symmetry_restored"] is False
