"""vle.py — Изотермический ЖПЭ-расчёт (равновесие пар-жидкость).

Порт MATLAB-функции ``matlab/functions/vle_flash.m``.
Уравнение Рэтчфорда-Райса + обобщённый закон Рауля с NRTL-поправкой.
Компонент C (ионная жидкость) считается нелетучим.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .nrtl import nrtl_gamma, R_GAS


@dataclass
class FlashResult:
    """Результат изотермического flash-расчёта."""
    V_frac: float       # мольная доля пара
    y: np.ndarray       # состав пара [yA, yB, 0]
    x: np.ndarray       # состав жидкости [xA, xB, xC]
    T_bub: float        # температура пузырькования [К]
    converged: bool


def psat(ant_row: np.ndarray, T: float) -> float:
    """Давление насыщенного пара по Антуану: ln(P[Па]) = A + B/(T+C)."""
    return float(np.exp(ant_row[0] + ant_row[1] / (T + ant_row[2])))


def bubble_T(
    z: np.ndarray,
    P: float,
    ant: np.ndarray,
    delta_g: np.ndarray,
    alpha: float,
    R: float = R_GAS,
    *,
    tol: float = 1e-9,
    max_iter: int = 500,
) -> float:
    """Температура пузырькования методом Ньютона-Рафсона.

    Условие: sum_i(K_i * z_i) = 1 при V=0, x=z.
    """
    z = np.asarray(z, dtype=float)
    T = 320.0
    for _ in range(max_iter):
        Ps = np.array([psat(ant[0], T), psat(ant[1], T)])
        g = nrtl_gamma(z, delta_g, alpha, T, R)
        f = g[0] * Ps[0] / P * z[0] + g[1] * Ps[1] / P * z[1] - 1.0

        dT = 0.1
        Ps2 = np.array([psat(ant[0], T + dT), psat(ant[1], T + dT)])
        g2 = nrtl_gamma(z, delta_g, alpha, T + dT, R)
        f2 = g2[0] * Ps2[0] / P * z[0] + g2[1] * Ps2[1] / P * z[1] - 1.0
        df = (f2 - f) / dT

        if abs(df) < 1e-15:
            break
        T_new = T - f / df
        if abs(T_new - T) < tol:
            T = T_new
            break
        T = T_new
    return T


def vle_flash(
    z: np.ndarray,
    P: float,
    T: float,
    ant: np.ndarray,
    delta_g: np.ndarray,
    alpha: float,
    R: float = R_GAS,
    *,
    tol: float = 1e-9,
    max_iter: int = 500,
) -> FlashResult:
    """Изотермический flash при заданных P, T.

    K_i = gamma_i * Psat_i / P (обобщённый закон Рауля).
    Уравнение Рэтчфорда-Райса решается бисекцией по V.
    """
    z = np.asarray(z, dtype=float).ravel()
    z = np.maximum(z, 1e-15)
    z = z / z.sum()

    Psat_vec = np.array([psat(ant[0], T), psat(ant[1], T), 0.0])

    # Начальные K (gamma при x≈z)
    x0 = z.copy()
    x0[2] = max(z[2], 1e-10)
    gamma0 = nrtl_gamma(x0, delta_g, alpha, T, R)
    K = gamma0 * Psat_vec / P
    K[2] = 0.0  # ИЖ не испаряется

    T_b = bubble_T(z, P, ant, delta_g, alpha, R, tol=tol, max_iter=max_iter)

    # Проверка однофазности (только жидкость)
    if np.sum(K * z) < 1.0:
        return FlashResult(0.0, np.zeros(3), z.copy(), T_b, True)

    V_lo, V_hi, V = 0.0, 1.0, 0.5
    converged = False
    x_it = z.copy()

    for _ in range(max_iter):
        denom = 1 + V * (K - 1)
        denom = np.where(denom < 1e-15, 1e-15, denom)
        x_it = z / denom
        x_it[2] = max(x_it[2], 0.0)
        if x_it.sum() > 0:
            x_it = x_it / x_it.sum()

        gamma_it = nrtl_gamma(x_it, delta_g, alpha, T, R)
        K = gamma_it * Psat_vec / P
        K[2] = 0.0

        RR = np.sum(z * (K - 1) / (1 + V * (K - 1)))
        if abs(RR) < tol:
            converged = True
            break

        if RR > 0:
            V_lo = V
        else:
            V_hi = V
        V = 0.5 * (V_lo + V_hi)

    x = np.maximum(x_it, 0)
    x = x / x.sum()
    y = K * x
    y[2] = 0.0
    if y.sum() > 0:
        y = y / y.sum()

    return FlashResult(min(max(V, 0.0), 1.0), y, x, T_b, converged)
