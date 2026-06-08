"""lle.py — Расчёт жидкость-жидкость равновесия (ЖЖЭ).

Порт MATLAB-функции ``matlab/functions/lle_solver.m``.
Метод последовательной подстановки с демпфированием.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .nrtl import nrtl_gamma, R_GAS


@dataclass
class LLEResult:
    """Результат расчёта ЖЖЭ."""
    xE: np.ndarray   # состав экстрактной фазы (богатой ИЖ)
    xR: np.ndarray   # состав рафинатной фазы (богатой MTBE)
    beta: float      # мольная доля экстрактной фазы
    iters: int       # число итераций
    flag: int        # 0: сошлось нетривиально; 1: не сошлось; 2: тривиально


def lle_solver(
    z: np.ndarray,
    delta_g: np.ndarray,
    alpha: float,
    T: float,
    R: float = R_GAS,
    *,
    tol: float = 1e-10,
    max_iter: int = 2000,
    damp: float = 0.2,
    xR0: np.ndarray | None = None,
    xE0: np.ndarray | None = None,
) -> LLEResult:
    """ЖЖЭ методом последовательной подстановки.

    Уравнение равновесия: gamma_i^R * x_i^R = gamma_i^E * x_i^E
    => K_i = x_i^R / x_i^E = gamma_i^E / gamma_i^R

    Демпфирование для устойчивости:
        xE_new = xE / sum(K .* xE)
        xE = damp*xE_new + (1-damp)*xE_old
    """
    z = np.asarray(z, dtype=float).ravel()
    z = z / z.sum()

    xR = np.array(xR0 if xR0 is not None else [0.97, 0.02, 0.01], dtype=float)
    xE = np.array(xE0 if xE0 is not None else [0.01, 0.02, 0.97], dtype=float)
    xR /= xR.sum()
    xE /= xE.sum()

    diff = 1.0
    it = 0
    flag = 1

    while diff > tol and it < max_iter:
        it += 1
        gR = nrtl_gamma(xR, delta_g, alpha, T, R)
        gE = nrtl_gamma(xE, delta_g, alpha, T, R)

        K = gE / gR  # коэффициент распределения x_R/x_E

        xR_old, xE_old = xR.copy(), xE.copy()

        xE_new = xE / np.sum(K * xE)
        xE_new = np.maximum(xE_new, 0)
        xE_new /= xE_new.sum()

        xE = damp * xE_new + (1 - damp) * xE_old
        xE = np.maximum(xE, 0)
        xE /= xE.sum()

        xR = K * xE
        xR = np.maximum(xR, 0)
        xR /= xR.sum()

        diff = np.max(np.abs(xR - xR_old)) + np.max(np.abs(xE - xE_old))

    beta = _phase_fraction(z, xR, xE)

    sep = np.sum(np.abs(xR - xE))
    if diff <= tol:
        flag = 0 if sep > 5e-3 else 2

    return LLEResult(xE=xE, xR=xR, beta=beta, iters=it, flag=flag)


def _phase_fraction(z: np.ndarray, xR: np.ndarray, xE: np.ndarray) -> float:
    """Доля экстрактной фазы: z = beta*xE + (1-beta)*xR."""
    d = xE - xR
    mask = np.abs(d) > 1e-6
    if np.any(mask):
        beta = float(np.mean((z[mask] - xR[mask]) / d[mask]))
        return min(max(beta, 0.0), 1.0)
    return 0.5
