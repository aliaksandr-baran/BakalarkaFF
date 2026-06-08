"""nrtl.py — Модель активностных коэффициентов NRTL.

Порт MATLAB-функции ``matlab/functions/nrtl_gamma.m``.
Система: MTBE (A) – Метанол (B) – [BMIM][HSO4] (C).
"""
from __future__ import annotations

import numpy as np

R_GAS = 8.314  # газовая постоянная [Дж/моль/К]


def nrtl_gamma(
    x: np.ndarray,
    delta_g: np.ndarray,
    alpha: float | np.ndarray,
    T: float,
    R: float = R_GAS,
) -> np.ndarray:
    """Коэффициенты активности по модели NRTL для n-компонентной смеси.

    Параметры
    ----------
    x : array_like [n]
        Мольные доли (нормируются автоматически).
    delta_g : ndarray [n, n]
        Матрица энергий взаимодействия [кДж/моль], diag = 0.
    alpha : float | ndarray
        Параметр нерандомности (скаляр или матрица n×n).
    T : float
        Температура [К].
    R : float
        Газовая постоянная [Дж/моль/К].

    Возвращает
    ----------
    gamma : ndarray [n]
        Коэффициенты активности.

    Уравнения NRTL (Chen, 1982; уравнения 14-17 диссертации):
        tau_ij = delta_g_ij * 1000 / (R*T)
        G_ij   = exp(-alpha * tau_ij)
        S_j    = sum_k(x_k * G_kj)
        ln_gamma_i = sum_k(x_k*tau_ki*G_ki)/S_i
                   + sum_j[x_j*G_ij/S_j * (tau_ij - sum_m(x_m*tau_mj*G_mj)/S_j)]
    """
    x = np.asarray(x, dtype=float).ravel()
    x = np.maximum(x, 1e-14)
    x = x / x.sum()
    n = x.size

    delta_g = np.asarray(delta_g, dtype=float)
    tau = (delta_g * 1000.0) / (R * T)       # [n, n]
    G = np.exp(-np.asarray(alpha) * tau)     # [n, n]

    # Нормировочные суммы S_j = sum_k(x_k * G_kj)
    S = x @ G                                 # [n]

    # numj = sum_m(x_m * tau_mj * G_mj) — числитель внутренней суммы
    numj = (x[:, None] * tau * G).sum(axis=0)   # [n]

    ln_gamma = np.zeros(n)
    for i in range(n):
        # Первое слагаемое: sum_k(x_k*tau_ki*G_ki) / S_i
        term1 = (x * tau[:, i] * G[:, i]).sum() / S[i]
        # Второе слагаемое: sum_j[x_j*G_ij/S_j * (tau_ij - numj/S_j)]
        term2 = np.sum(
            (x * G[i, :] / S) * (tau[i, :] - numj / S)
        )
        ln_gamma[i] = term1 + term2

    return np.exp(ln_gamma)
