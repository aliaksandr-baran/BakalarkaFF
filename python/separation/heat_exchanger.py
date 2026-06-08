"""heat_exchanger.py — Тепловой расчёт теплообменника (метод LMTD).

Порт MATLAB-функции ``matlab/functions/heat_exchanger.m``.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class HXResult:
    Q: float
    LMTD_cf: float
    F_corr: float
    LMTD: float
    U_design: float
    A_req: float
    n_tubes: int
    n_plates: int
    effectiveness: float
    hx_type: str


def _select_U(hx_type: str, Th_in: float) -> float:
    """Типовые значения U [Вт/м²/К] (Kern, 1950)."""
    if hx_type == "plate":
        return 3000.0
    return 1200.0 if abs(Th_in - 145) < 30 else 800.0


def _lmtd_correction(Th_in, Th_out, Tc_in, Tc_out, passes) -> float:
    """Поправочный коэффициент F (TEMA). F≈1 для противотока/фазового перехода."""
    if passes == 1:
        return 1.0
    Rr = (Th_in - Th_out) / max(Tc_out - Tc_in, 1e-6)
    P = (Tc_out - Tc_in) / max(Th_in - Tc_in, 1e-6)
    S = np.sqrt(Rr ** 2 + 1)
    arg = (2 / P - 1 - Rr + S) / (2 / P - 1 - Rr - S)
    if arg <= 0 or arg == 1 or abs(Rr - 1) < 1e-4:
        return 1.0
    F = S * np.log((1 - P * Rr) / max(1 - P, 1e-9)) / \
        ((Rr - 1) * np.log(max(arg, 1e-9)))
    return min(max(F, 0.5), 1.0)


def heat_exchanger(
    Th_in: float, Th_out: float,
    Tc_in: float, Tc_out: float,
    Q_duty: float,
    *,
    hx_type: str = "shell_tube",
    U: float | None = None,
    passes: int = 1,
    d_tube: float = 0.025,
    L_tube: float = 3.0,
    fouling_R: float = 1e-4,
    name: str = "Теплообменник",
    verbose: bool = True,
) -> HXResult:
    """Тепловой расчёт по LMTD. Температуры в °C, Q в кВт."""
    if U is None:
        U = _select_U(hx_type, Th_in)

    dT1, dT2 = Th_in - Tc_out, Th_out - Tc_in
    if abs(dT1 - dT2) < 1e-4:
        LMTD_cf = dT1
    elif dT1 <= 0 or dT2 <= 0:
        LMTD_cf = float("nan")
    else:
        LMTD_cf = (dT1 - dT2) / np.log(dT1 / dT2)

    F = _lmtd_correction(Th_in, Th_out, Tc_in, Tc_out, passes)
    LMTD = LMTD_cf * F
    U_eff = 1.0 / (1.0 / U + fouling_R)

    A_req = (Q_duty * 1000) / (U_eff * LMTD) if (LMTD and LMTD > 0) else float("nan")

    A_tube = np.pi * d_tube * L_tube
    n_tubes = int(np.ceil(A_req / A_tube)) if np.isfinite(A_req) else 0
    A_plate = 0.5 * 1.5
    n_plates = int(np.ceil(A_req / A_plate)) + 1 if np.isfinite(A_req) else 0

    if abs(Th_in - Th_out) < 0.5:
        eff = 1.0
    else:
        eff = (Th_in - Th_out) / max(Th_in - Tc_in, 1e-3)

    if verbose:
        extra = f"n_труб={n_tubes}" if hx_type == "shell_tube" else f"n_пластин={n_plates}"
        print(f"  {name:<28}  Q={Q_duty:7.1f} кВт  LMTD={LMTD:6.2f} K  "
              f"U={U_eff:5.0f}  A={A_req:6.2f} м²  {extra}")

    return HXResult(
        Q=Q_duty, LMTD_cf=LMTD_cf, F_corr=F, LMTD=LMTD, U_design=U_eff,
        A_req=A_req, n_tubes=n_tubes, n_plates=n_plates,
        effectiveness=eff, hx_type=hx_type,
    )
