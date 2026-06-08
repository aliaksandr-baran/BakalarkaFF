"""Базовые тесты пакета separation (проверка паритета с MATLAB)."""
import numpy as np
import pytest

import separation as sep
from separation import thermo


def test_nrtl_gamma_normalization():
    """gamma > 0 и для чистого компонента gamma -> 1."""
    g = sep.nrtl_gamma([1.0, 0.0, 0.0], thermo.DELTA_G, thermo.ALPHA, 318.15)
    assert np.all(g > 0)
    assert g[0] == pytest.approx(1.0, abs=1e-6)


def test_nrtl_gamma_infinite_dilution():
    """В бесконечном разбавлении gamma конечна и положительна."""
    g = sep.nrtl_gamma([0.001, 0.5, 0.499], thermo.DELTA_G, thermo.ALPHA, 318.15)
    assert np.all(np.isfinite(g))
    assert np.all(g > 0)


def test_lle_two_phase_split():
    """В двухфазной области фазы должны различаться (flag=0)."""
    res = sep.lle_solver([0.3, 0.3, 0.4], thermo.DELTA_G, thermo.ALPHA, 318.15)
    assert res.flag == 0
    assert np.sum(np.abs(res.xR - res.xE)) > 5e-3
    assert res.xR[0] > res.xE[0]   # рафинат богаче MTBE


def test_vle_flash_nonvolatile_il():
    """ИЖ (компонент C) не должна попадать в пар."""
    fr = sep.vle_flash([0.18, 0.82, 0.0], 50e3, 318.0,
                       thermo.ANT, thermo.DELTA_G, thermo.ALPHA)
    assert fr.y[2] == 0.0


def test_psat_monotonic():
    """Давление насыщения растёт с температурой."""
    p1 = sep.psat(thermo.ANT[1], 300.0)
    p2 = sep.psat(thermo.ANT[1], 330.0)
    assert p2 > p1


def test_pump_positive_power():
    """При положительном перепаде давления мощность положительна."""
    r = sep.pump_sizing(100, [0.5, 0.5, 0.0], 101325, 301325,
                        dz=2, verbose=False)
    assert r.P_shaft > 0
    assert r.Q_vol > 0


def test_heat_exchanger_area():
    """Площадь теплообмена положительна для корректных температур."""
    r = sep.heat_exchanger(145, 145, 25, 45, 437.02, U=750, verbose=False)
    assert r.A_req > 0
    assert r.n_tubes > 0


def test_economics_capex_indexing():
    """CEPCI-индексация увеличивает стоимость при росте индекса."""
    eq = [sep.Equipment("test", cost_base=1000, CEPCI_base=567.5)]
    util = sep.UtilityRates(12, 0.693, 0.172, 2.8, 73.67, 4.96)
    res = sep.process_economics(eq, util, 7920, CEPCI_now=811.5,
                                IL_cost=0, verbose=False)
    assert res.equip_costs[0] == pytest.approx(1000 * 811.5 / 567.5)
    assert res.TCI > res.CAPEX_equip


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-v"]))
