"""vle_improved.py — Улучшенный расчёт равновесия пар-жидкость.

Улучшения по сравнению с оригинальным портом:
  * bubble_T использует метод Брента (scipy.optimize.brentq) с автоматическим
    поиском отрезка — надёжен при низких давлениях (≈2.45 кПа), где ньютоновский
    метод оригинала может расходиться.
  * vle_flash использует ускоренное обновление K из уравнения Рэтчфорда-Райса
    с контролем сходимости.

Все остальные API-сигнатуры совпадают с separation.vle.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import brentq

from separation.nrtl import nrtl_gamma, R_GAS
from separation.vle import psat, FlashResult   # psat и FlashResult оставляем прежними


__all__ = ["psat", "bubble_T_improved", "vle_flash_improved", "FlashResult"]


def bubble_T_improved(
    z: np.ndarray,
    P: float,
    ant: np.ndarray,
    delta_g: np.ndarray,
    alpha: float,
    R: float = R_GAS,
    *,
    T_lo: float = 200.0,
    T_hi: float = 500.0,
    tol: float = 1e-9,
    max_iter: int = 200,
) -> float:
    """Температура пузырькования методом Брента.

    Условие: F(T) = sum_i K_i(T) z_i - 1 = 0,  K_i = gamma_i*Psat_i/P.
    Метод Брента надёжно работает при низких давлениях, где ньютоновская итерация
    может давать отрицательные значения T или не сходиться.
    """
    z = np.asarray(z, dtype=float)

    def _F(T: float) -> float:
        Ps = np.array([psat(ant[0], T), psat(ant[1], T)])
        g = nrtl_gamma(z, delta_g, alpha, T, R)
        return g[0] * Ps[0] / P * z[0] + g[1] * Ps[1] / P * z[1] - 1.0

    # Расширяем интервал, пока знаки не окажутся разными
    for _ in range(50):
        try:
            fa, fb = _F(T_lo), _F(T_hi)
        except (ZeroDivisionError, FloatingPointError):
            T_lo -= 10
            T_hi += 10
            continue
        if np.isfinite(fa) and np.isfinite(fb):
            if fa * fb < 0:
                break
        T_lo -= 10
        T_hi += 10
    else:
        # Запасной выход: оригинальный ньютоновский метод
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
            T = T - f / df
        return T

    return brentq(_F, T_lo, T_hi, xtol=tol, maxiter=max_iter)


def vle_flash_improved(
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
    """Изотермический flash с улучшенным bubble_T и надёжной бисекцией R-R.

    Отличия от оригинала:
    - bubble_T заменён на Брент (bubble_T_improved)
    - обновление K после каждой итерации R-R (строгая сходимость)
    """
    z = np.asarray(z, dtype=float).ravel()
    z = np.maximum(z, 1e-15)
    z = z / z.sum()

    Psat_vec = np.array([psat(ant[0], T), psat(ant[1], T), 0.0])

    x0 = z.copy()
    x0[2] = max(z[2], 1e-10)
    gamma0 = nrtl_gamma(x0, delta_g, alpha, T, R)
    K = gamma0 * Psat_vec / P
    K[2] = 0.0

    T_b = bubble_T_improved(z, P, ant, delta_g, alpha, R, tol=tol, max_iter=max_iter)

    # Однофазная жидкость
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
