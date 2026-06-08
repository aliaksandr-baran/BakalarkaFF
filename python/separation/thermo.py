"""thermo.py — Физические свойства и термодинамические корреляции.

Теплоёмкости, теплоты испарения, плотности компонентов.
Общие параметры для функций теплового расчёта.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

# Молярные массы [кг/кмоль]
M_A = 88.15    # MTBE
M_B = 32.04    # Метанол
M_C = 236.29   # [BMIM][HSO4]
M = np.array([M_A, M_B, M_C])

# --- Параметры NRTL (T = 318.15 K) ---
ALPHA = 0.4
DELTA_G = np.array([
    [0.0000, 5.5737, 19.2694],
    [7.4991, 0.0000, 10.6456],
    [16.2120, 0.5632, 0.0000],
])  # [кДж/моль]

# --- Параметры Антуана: ln(Psat[Па]) = A + B/(T+C), T[К] ---
ANT = np.array([
    [20.84054, -2624.525, -46.15171],   # A = MTBE
    [23.5347, -3661.468, -32.77001],    # B = MeOH
])

# --- Коэффициенты теплоёмкости ---
CP_A_KGK = 2.127                                         # MTBE [кДж/кг/К]
CP_B_COEFF = [0.8382, -0.003231, 8.296e-6, -1.689e-10]   # MeOH DIPPR [ккал/кг/К]
CP_IL_COEFF = [753.28, -3.4195, 0.0063]                  # ИЖ [кДж/кмоль/К]


@dataclass
class VapParams:
    """Параметры уравнения Ватсона-Масса для теплоты испарения."""
    A: float
    Tc: float
    beta: float
    alpha: float = 0.0   # для MeOH (показатель зависит от T)


VAP_A = VapParams(A=46.23, Tc=497.1, beta=0.2893)               # MTBE
VAP_B = VapParams(A=45.3, Tc=512.6, beta=0.4241, alpha=-0.31)   # MeOH


def cp_mtbe(T: float) -> float:
    """Теплоёмкость MTBE [кДж/кмоль/К] (постоянная)."""
    return CP_A_KGK * M_A


def cp_methanol(T: float) -> float:
    """Теплоёмкость метанола [кДж/кмоль/К] (DIPPR), T[К]."""
    a, b, c, d = CP_B_COEFF
    return (a + b * T + c * T ** 2 + d * T ** 3) * M_B * 4.184


def cp_il(T: float) -> float:
    """Теплоёмкость ИЖ [кДж/кмоль/К], t = T-273.15."""
    a, b, c = CP_IL_COEFF
    t = T - 273.15
    return a + b * t + c * t ** 2


def cp_mix(x: np.ndarray, T: float) -> float:
    """Среднесмесевая теплоёмкость жидкости [кДж/кмоль/К]."""
    x = np.asarray(x)
    return x[0] * cp_mtbe(T) + x[1] * cp_methanol(T) + x[2] * cp_il(T)


def hvap_mtbe(T: float, vap: VapParams = VAP_A) -> float:
    """Теплота испарения MTBE [кДж/моль] (Ватсон-Масс)."""
    base = max((vap.Tc - T) / (vap.Tc - 298.15), 0.0)
    return vap.A * base ** vap.beta


def hvap_methanol(T: float, vap: VapParams = VAP_B) -> float:
    """Теплота испарения метанола [кДж/моль] (Ватсон-Масс)."""
    base = max((vap.Tc - T) / (vap.Tc - 298.15), 0.0)
    exponent = vap.alpha + vap.beta * (T / vap.Tc)
    return vap.A * base ** exponent
