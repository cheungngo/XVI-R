import json
from pathlib import Path
R = Path(__file__).resolve().parents[2]


def j(n):
    return json.loads((R / "certificates" / f"{n}.json").read_text())


def test_gram():
    g = j("R08_gram")
    assert g["G_nonzero"] and g["G_bose_symmetric"] and g["G_sigma_image_zero"] and g["G_in_4D_basis"]
    assert g["H_gram_reduced"] == g["H_formal"] == 6 and g["G_blind"] and not g["G_trivial"]
    assert [k for k, v in g["slot_census"].items() if v["gram_possible_4D"]] == ["H"]


def test_parity_and_fourpoint_honest():
    assert j("R08_parity")["odd_cohomology_dimension"].startswith("Missing")
    assert j("R08_four_point")["outcome"].startswith("OPEN")
    assert j("R08_P6_status")["bulk_symmetry_restored"] is False
