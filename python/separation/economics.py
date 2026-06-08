"""economics.py — CAPEX / OPEX расчёт с индексированием CEPCI.

Порт MATLAB-функции ``matlab/functions/process_economics.m``.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Sequence

import numpy as np


@dataclass
class Equipment:
    name: str
    cost_base: float           # базовая стоимость [€]
    CEPCI_base: float = 567.5  # 2017
    eq_type: str = "other"


@dataclass
class UtilityRates:
    steam_price: float    # €/т
    cooling_price: float  # €/т
    elec_price: float     # €/кВт·ч
    steam_kg_h: float     # т/ч
    cooling_kg_h: float   # т/ч
    elec_kW: float        # кВт


@dataclass
class EconomicsResult:
    equip_costs: np.ndarray
    CAPEX_equip: float
    CAPEX_install: float
    TCI: float
    OPEX: float
    OPEX_breakdown: dict


def process_economics(
    equipment: Sequence[Equipment],
    util: UtilityRates,
    annual_hours: float,
    *,
    CEPCI_now: float = 811.5,
    f_install: float = 3.1,
    f_TCI: float = 1.18,
    f_labor: float = 0.02,
    f_maint: float = 0.01,
    IL_cost: float = 0.0,
    verbose: bool = True,
) -> EconomicsResult:
    """CAPEX (индекс CEPCI + Lang-фактор) и OPEX (энергоносители + труд + ТО)."""
    costs = np.array([
        eq.cost_base * (CEPCI_now / eq.CEPCI_base) for eq in equipment
    ])

    CAPEX_equip = costs.sum() + IL_cost
    CAPEX_install = CAPEX_equip * f_install
    TCI = CAPEX_install * f_TCI

    cost_steam = util.steam_price * util.steam_kg_h * annual_hours / 1000
    cost_cooling = util.cooling_price * util.cooling_kg_h * annual_hours / 1000
    cost_elec = util.elec_price * util.elec_kW * annual_hours
    cost_labor = TCI * f_labor
    cost_maint = TCI * f_maint
    OPEX = cost_steam + cost_cooling + cost_elec + cost_labor + cost_maint

    breakdown = {
        "steam": cost_steam, "cooling": cost_cooling, "elec": cost_elec,
        "labor": cost_labor, "maint": cost_maint,
    }

    if verbose:
        print("\n--- CAPEX: Стоимость оборудования ---")
        print(f"{'Оборудование':<22} {'C_base [€]':>12} {'CEPCI':>8} {'C_2025 [€]':>12}")
        print("-" * 58)
        for eq, c in zip(equipment, costs):
            print(f"{eq.name:<22} {eq.cost_base:>12.0f} {eq.CEPCI_base:>8.1f} {c:>12.0f}")
        print("-" * 58)
        print(f"  Итого оборудование:  {CAPEX_equip:>12.0f} €")
        print(f"  С монтажом (f={f_install}): {CAPEX_install:>12.0f} €")
        print(f"  TCI:                 {TCI:>12.0f} €")
        print("\n--- OPEX: Операционные расходы [€/год] ---")
        print(f"  Греющий пар:         {cost_steam:>12.0f}")
        print(f"  Охлаждающая вода:    {cost_cooling:>12.0f}")
        print(f"  Электроэнергия:      {cost_elec:>12.0f}")
        print(f"  Труд персонала:      {cost_labor:>12.0f}")
        print(f"  Техническое обслуж.: {cost_maint:>12.0f}")
        print(f"  ИТОГО OPEX:          {OPEX:>12.0f}  €/год")

    return EconomicsResult(
        equip_costs=costs, CAPEX_equip=CAPEX_equip,
        CAPEX_install=CAPEX_install, TCI=TCI, OPEX=OPEX,
        OPEX_breakdown=breakdown,
    )
