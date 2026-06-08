"""evaporator.py — Расчёт испарителя (однократное испарение + тепловой баланс).

Порт MATLAB-функции ``matlab/functions/evaporator_calc.m``.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import integrate

from . import thermo
from .nrtl import R_GAS
from .vle import vle_flash, bubble_T


@dataclass
class EvaporatorResult:
    vapor_n: float
    vapor_x: np.ndarray
    liquid_n: float
    liquid_x: np.ndarray
    T_op: float
    V_frac: float
    Q_duty: float        # кВт
    Q_sensible: float
    Q_latent: float
    A_evap: float        # м²
    LMTD: float


def _lmtd(Th_in, Th_out, Tc_in, Tc_out):
    dT1, dT2 = Th_in - Tc_out, Th_out - Tc_in
    if abs(dT1 - dT2) < 1e-4:
        return dT1
    if dT1 <= 0 or dT2 <= 0:
        return max(abs(dT1), abs(dT2))
    return (dT1 - dT2) / np.log(dT1 / dT2)


def evaporator_calc(
    n: float,
    x: np.ndarray,
    T_in: float,
    P_op: float,
    ant: np.ndarray,
    delta_g: np.ndarray,
    alpha: float,
    R: float = R_GAS,
    *,
    T_steam: float = 418.15,
    U_evap: float = 600.0,
    V_spec: float | None = None,
    name: str = "Испаритель",
    verbose: bool = True,
) -> EvaporatorResult:
    """Однократное равновесное испарение при давлении P_op.

    Тепловая нагрузка = теплота нагрева жидкости (Q_sensible)
    + теплота испарения пара (Q_latent).
    """
    x = np.asarray(x, dtype=float)
    x = np.maximum(x, 0)
    x = x / x.sum()

    T_op = bubble_T(x, P_op, ant, delta_g, alpha, R)
    fr = vle_flash(x, P_op, T_op, ant, delta_g, alpha, R)
    V_frac = fr.V_frac if V_spec is None else V_spec
    V_frac = min(max(V_frac, 0.0), 1.0)
    y_out, x_out = fr.y, fr.x

    n_vap = n * V_frac
    n_liq = n * (1 - V_frac)

    # Теплота нагрева жидкости [кДж/кмоль] -> кВт
    if abs(T_op - T_in) > 0.1:
        dH_sens, _ = integrate.quad(lambda T: thermo.cp_mix(x, T), T_in, T_op)
    else:
        dH_sens = 0.0
    Q_sensible = n * dH_sens / 3600.0

    # Теплота испарения [кДж/моль] -> кДж/кмоль
    Hvap_A = thermo.hvap_mtbe(T_op)
    Hvap_B = thermo.hvap_methanol(T_op)
    if y_out[:2].sum() > 1e-9:
        Hvap_mix = (y_out[0] * Hvap_A + y_out[1] * Hvap_B) * 1000
    else:
        Hvap_mix = 0.0
    Q_latent = n_vap * Hvap_mix / 3600.0

    Q_total = Q_sensible + Q_latent

    LMTD = _lmtd(T_steam - 273.15, T_steam - 273.15, T_in - 273.15, T_op - 273.15)
    A_evap = Q_total * 1000 / (U_evap * LMTD) if LMTD > 0 else float("nan")

    if verbose:
        print(f"  {name:<20}  P={P_op/1e3:6.2f} кПа  "
              f"T_кип={T_op-273.15:6.2f}°C  V={V_frac:.3f}")
        print(f"    Q_сенс={Q_sensible:7.2f} кВт  Q_лат={Q_latent:7.2f} кВт  "
              f"Q_total={Q_total:7.2f} кВт")
        print(f"    Площадь: {A_evap:.2f} м²   LMTD: {LMTD:.2f} К")
        print(f"    Пар:  n={n_vap:.2f}  xA={y_out[0]:.4f} xB={y_out[1]:.4f} xC={y_out[2]:.4f}")
        print(f"    Жидк: n={n_liq:.2f}  xA={x_out[0]:.4f} xB={x_out[1]:.4f} xC={x_out[2]:.4f}")

    return EvaporatorResult(
        vapor_n=n_vap, vapor_x=y_out, liquid_n=n_liq, liquid_x=x_out,
        T_op=T_op, V_frac=V_frac, Q_duty=Q_total,
        Q_sensible=Q_sensible, Q_latent=Q_latent, A_evap=A_evap, LMTD=LMTD,
    )
