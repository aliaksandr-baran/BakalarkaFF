"""pump.py — Расчёт насоса: гидравлика, мощность, подбор типа.

Порт MATLAB-функции ``matlab/functions/pump_sizing.m``.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import thermo

# Плотности компонентов [кг/м³] при ~45°C
RHO = np.array([740.0, 775.0, 1275.0])   # MTBE, MeOH, [BMIM][HSO4]
G = 9.81


@dataclass
class PumpResult:
    Q_vol: float          # м³/ч
    rho_mix: float        # кг/м³
    dP: float             # Па
    H_head: float         # м
    P_hydraulic: float    # Вт
    P_shaft: float        # Вт
    P_motor: float        # кВт
    NPSH_req: float       # м
    pump_type: str


def _pump_type(Q_vol: float, dP: float) -> str:
    dP_bar = dP / 1e5
    if Q_vol < 5 and dP_bar > 10:
        return "плунжерный"
    if Q_vol < 20:
        return "центробежный (малый)"
    if Q_vol < 200:
        return "центробежный"
    return "центробежный (крупный)"


def pump_sizing(
    n: float,
    x: np.ndarray,
    P_in: float,
    P_out: float,
    *,
    eta_pump: float = 0.70,
    eta_motor: float = 0.95,
    dz: float = 1.0,
    rho: np.ndarray = RHO,
    name: str = "Насос",
    verbose: bool = True,
) -> PumpResult:
    """Гидравлический расчёт насоса. Плотность смеси по правилу объёмов."""
    x = np.asarray(x, dtype=float)
    x = np.maximum(x, 0)
    x = x / x.sum()

    M_mix = float(x @ thermo.M)
    w = x * thermo.M / M_mix                    # массовые доли
    rho_mix = 1.0 / np.sum(w / rho)

    G_mass = n * M_mix                          # кг/ч
    Q_vol = G_mass / rho_mix                    # м³/ч
    Q_vol_s = Q_vol / 3600.0

    dP = P_out - P_in
    H_head = dP / (rho_mix * G) + dz

    P_hydr = Q_vol_s * dP
    P_shaft = P_hydr / eta_pump
    P_motor = P_shaft / eta_motor / 1000.0

    NPSH_req = 0.3 + 0.03 * np.sqrt(Q_vol)
    ptype = _pump_type(Q_vol, dP)

    if verbose:
        print(f"  {name:<16}  Q={Q_vol:6.2f} м³/ч  dP={dP/1e3:7.1f} кПа  "
              f"H={H_head:6.1f} м  P_вал={P_shaft:7.1f} Вт  "
              f"P_дв={P_motor:5.2f} кВт  [{ptype}]")

    return PumpResult(
        Q_vol=Q_vol, rho_mix=rho_mix, dP=dP, H_head=H_head,
        P_hydraulic=P_hydr, P_shaft=P_shaft, P_motor=P_motor,
        NPSH_req=NPSH_req, pump_type=ptype,
    )
