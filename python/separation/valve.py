"""valve.py — Изоэнтальпийное дросселирование (клапан).

Порт MATLAB-функции ``matlab/functions/throttle_valve.m``.
H_in = H_out: при понижении давления возможно частичное испарение.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import thermo
from .nrtl import R_GAS
from .vle import vle_flash, bubble_T


@dataclass
class ValveResult:
    T_out: float
    P_out: float
    phase: str
    V_frac: float
    T_bub: float
    dT: float
    joule_thomson: float
    vapor_n: float
    vapor_x: np.ndarray
    liquid_n: float
    liquid_x: np.ndarray


def throttle_valve(
    n: float,
    x: np.ndarray,
    T_in: float,
    P_in: float,
    P_out: float,
    ant: np.ndarray,
    delta_g: np.ndarray,
    alpha: float,
    R: float = R_GAS,
    *,
    name: str = "Клапан",
    verbose: bool = True,
) -> ValveResult:
    """Изоэнтальпийное дросселирование с проверкой фазового состояния.

    Если T_in <= T_пузырьк(P_out): поток остаётся жидкостью.
    Иначе решается изоэнтальпийный баланс для двухфазной смеси:
        cp_mix*(T_in - T_out) = V_frac * Hvap_mix(T_out)
    """
    x = np.asarray(x, dtype=float)
    x = np.maximum(x, 0)
    x = x / x.sum()

    def hvap_mix(T, y):
        return y[0] * thermo.hvap_mtbe(T) * 1000 + y[1] * thermo.hvap_methanol(T) * 1000

    T_bub = bubble_T(x, P_out, ant, delta_g, alpha, R)

    if T_in <= T_bub + 0.5:
        T_out, V_frac, phase = T_in, 0.0, "liquid"
        y_out, x_out = np.zeros(3), x.copy()
    else:
        T_out = T_bub
        y_out, x_out, V_frac = np.zeros(3), x.copy(), 0.0
        for _ in range(100):
            fr = vle_flash(x, P_out, T_out, ant, delta_g, alpha, R)
            if fr.V_frac <= 0:
                T_out, V_frac = T_bub, 0.0
                y_out, x_out = np.zeros(3), x.copy()
                break
            Hv = hvap_mix(T_out, fr.y)
            cp_val = 0.5 * (thermo.cp_mix(x, T_in) + thermo.cp_mix(x, T_out))
            V_new = min(max(cp_val * (T_in - T_out) / max(Hv, 1), 0.0), 1.0)
            if abs(V_new - fr.V_frac) < 1e-5:
                V_frac, y_out, x_out = fr.V_frac, fr.y, fr.x
                break
            T_out_new = T_in - V_new * Hv / max(cp_val, 1)
            T_out_new = max(min(T_out_new, T_in), T_bub - 5)
            T_out = 0.5 * T_out_new + 0.5 * T_out
        T_out = max(T_out, T_bub)
        phase = "vapor" if V_frac > 0.99 else "two-phase"

    dP = P_in - P_out
    mu_JT = (T_in - T_out) / dP if abs(dP) > 1e3 else 0.0

    n_vap = n * V_frac
    n_liq = n * (1 - V_frac)

    if verbose:
        print(f"  {name:<20}  P:{P_in/1e3:6.1f}→{P_out/1e3:6.1f} кПа  "
              f"T_in={T_in-273.15:6.2f}°C  T_out={T_out-273.15:6.2f}°C  "
              f"ΔT={T_in-T_out:.3f} К")
        print(f"    Фаза: {phase:<10}  V_frac={V_frac:.4f}  "
              f"n_пар={n_vap:.2f}  n_жидк={n_liq:.2f} кмоль/ч   μ_ДТ={mu_JT:.3e} К/Па")

    return ValveResult(
        T_out=T_out, P_out=P_out, phase=phase, V_frac=V_frac, T_bub=T_bub,
        dT=T_in - T_out, joule_thomson=mu_JT,
        vapor_n=n_vap, vapor_x=y_out, liquid_n=n_liq, liquid_x=x_out,
    )
