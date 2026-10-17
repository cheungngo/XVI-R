import json
from pathlib import Path
import sympy as sp
R = Path(__file__).resolve().parents[2]


def j(n):
    return json.loads((R / "certificates" / f"{n}.json").read_text())


def test_hat_rule_unit():
    import sys
    for p in ("r07/code", "r06/code", "r05/code", "r03/code"):
        sys.path.insert(0, str(R / p))
    import r07_hatgen as HG
    T, H = sp.Symbol, sp.Symbol
    e = T("T_KU_c") * (H("H_EU_x") + T("T_EU_x"))
    assert HG.contract_hatgen(e, [("KU", "EU")]) == T("T_c_x")
    e2 = H("H_KU_y") * (H("H_EU_x") + T("T_EU_x"))
    assert HG.contract_hatgen(e2, [("KU", "EU")]) == H("H_x_y")
    e3 = (H("H_KU_y") + T("T_KU_y")) * (H("H_EU_x") + T("T_EU_x"))
    assert HG.contract_hatgen(e3, [("KU", "EU")]) == H("H_x_y") + T("T_x_y")


def test_sanity_and_classification():
    c = j("R07_classification")
    assert c["R05rule_sanity"]["dim_H"] == 18 and c["R05rule_sanity"]["dim_Nbar_g"] == 133
    h = c["hatgen"]
    assert h["sigma_inside_Nbar"] and h["images_inside_Nbar"]
    assert c["H_dim"] == h["dim_H"] == j("R07_H_basis")["dim"]


def test_front1():
    f = j("R07_front1_reduced_hhc")
    assert f["R05rule"]["full_minus_reduced_equals_Dgf_plus_sigmaE"] and f["R05rule"]["tree"]["A_plus_B_EH_zero"]
    assert f["outcome"] == "CLOSED" and f["hatgen"]["Rbar_4D_zero"]


def test_p22():
    p = j("R07_p22")
    for g in ("D01", "D03", "D11"):
        assert p["identification"][g]["same_ratio"] and p["identification"][g]["p22_CT_equals_minus_pole"]
    assert p["completion"]["Xi_4D_zero"] and p["outcome"].startswith("CLOSED")


def test_status():
    s = j("R07_P6_status")
    assert s["P6_outcome"] == "Partial" and s["bulk_symmetry_restored"] is False
