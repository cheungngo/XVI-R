import json
from pathlib import Path
R = Path(__file__).resolve().parents[2]


def j(n):
    return json.loads((R / "certificates" / f"{n}.json").read_text())


def test_A_changed_evanescent_and_fit_solvable():
    a = j("R06_loop_remainders")
    assert a["loop_remainders_unchanged"] is False and a["difference_purely_evanescent"] is True
    q = a["refit"]["Q_full_hhc_hat"]
    assert q["solvable_with_ghost_eq"] and q["particular_CT_all_evanescent"] and q["verified_rows_zero"]
    assert a["earlier_claims"]["Q07_all_four_rows_close_in_4D"].startswith("UNAFFECTED")
    assert a["earlier_claims"]["ghost_equation_compatibility"].startswith("UNAFFECTED")


def test_reduced_homogeneous_matches_R05():
    a = j("R06_loop_remainders")["refit"]["R_reduced_hhc_hat"]
    assert a["dim_N_with_ghost_eq"] == j("R06_classification")["Nbar_g"] == 133


def test_B_completeness():
    c = j("R06_completeness")
    assert all(v["partition_ok"] and v["frozen_equals_independent_ev"] for v in c["counts"].values())
    assert c["leak_outside_ev_plus_4d"] == 0 and c["pure_ev_images_inside_Nbar"] is True
    assert c["enlarged"]["dim_H"] == 18 and c["enlarged"]["new_trivial_directions"] == 0
    assert c["transport_R05_18"]["stable"] is True
    assert c["R05_reproduced"]["dim_H"] == 18


def test_status():
    s = j("R06_P6_status")
    assert s["P6_outcome"] == "Partial" and s["bulk_symmetry_restored"] is False
    assert j("R06_classification")["physical_certified"] == 0
