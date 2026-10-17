import json
from pathlib import Path
R = Path(__file__).resolve().parents[2]


def j(n):
    return json.loads((R / "certificates" / f"{n}.json").read_text())


def test_gap_closed_for_F_and_KKc():
    r = j("R03_reduced")
    assert r["F_images_failing_reduced"] == 0 and set(r["KKc"]["blocks"]) == {"pass"}


def test_diagnosis_localised():
    lg = j("R03_gap_diagnosis")["leg_identities_hhc_block"]
    assert lg["I1_EH_leg1"] == 0 and lg["I2_legs23_vs_A_K"] == 0 and lg["I3_E_legs23_vs_A_C"] == 0
    assert lg["I3_gf1_plus_E1"] == lg["total"]


def test_classification_adds_up():
    t = j("R03_classification")["transport_R02_22"]
    assert sum(t[k] for k in ("fail_reduced_constraints", "trivial_or_blind_in_reduced", "nontrivial_classes_in_H_bar")) == 22


def test_no_physical():
    assert j("R03_classification")["physical_certified"] == 0
