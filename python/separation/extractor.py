"""extractor.py — Противоточный жидкость-жидкость экстрактор.

Порт MATLAB-функции ``matlab/functions/extractor_column.m``.
Материальный баланс + масса-перенос + габариты колонны.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .nrtl import nrtl_gamma, R_GAS


@dataclass
class Stream:
    """Технологический поток."""
    n: float                # молярный поток [кмоль/ч]
    x: np.ndarray           # мольные доли [xA, xB, xC]
    T: float = 318.15       # температура [К]
    name: str = ""

    def __post_init__(self):
        self.x = np.asarray(self.x, dtype=float)


@dataclass
class ExtractorResult:
    raffinate: Stream
    extract: Stream
    stage_xR: np.ndarray
    stage_xE: np.ndarray
    stage_nR: np.ndarray
    stage_nE: np.ndarray
    K_profile: np.ndarray
    N_stages_req: int
    E_stage: float
    col_height: float
    col_diameter: float
    flooding_vel: float
    oper_vel: float
    K_La: float
    NTU: float
    HTU: float
    HETS: float
    iters: int
    error: float


# Молярные массы [кг/кмоль]
_M = np.array([88.15, 32.04, 236.29])


def extractor_column(
    feed: Stream,
    solvent: Stream,
    N_stages: int,
    delta_g: np.ndarray,
    alpha: float,
    T: float,
    R: float = R_GAS,
    *,
    tol: float = 1e-8,
    max_iter: int = 3000,
    damp: float = 0.15,
    E_stage: float = 0.85,
    HETS: float = 0.15,
    rho_R: float = 750.0,
    rho_E: float = 1230.0,
    sigma: float = 15.0,
    verbose: bool = True,
) -> ExtractorResult:
    """Расчёт противоточного ЖЖ-экстрактора с эффективностью Мёрфри.

    Нумерация: строка 0 = питание / выход экстракта,
               строка N = растворитель / выход рафината.
    """
    N = N_stages
    n_F, x_F = feed.n, feed.x
    n_S, x_S = solvent.n, solvent.x

    xR = np.zeros((N + 1, 3))
    xE = np.zeros((N + 1, 3))
    nR = np.zeros(N + 1)
    nE = np.zeros(N + 1)

    xR[0] = x_F
    nR[0] = n_F
    xE[N] = x_S
    nE[N] = n_S

    # Линейная интерполяция начального профиля
    for i in range(1, N + 1):
        t = i / N
        row = (1 - t) * x_F + t * np.array([0.99, 0.005, 0.005])
        xR[i] = row / row.sum()
    for i in range(N):
        t = (N - i) / N
        row = (1 - t) * x_S + t * np.array([0.04, 0.55, 0.41])
        xE[i] = row / row.sum()

    n_total = n_F + n_S
    x_total = (n_F * x_F + n_S * x_S) / n_total
    nE[0] = min(n_total * x_total[2] / max(xE[0, 2], 1e-6), n_total * 0.9)
    nR[N] = n_total - nE[0]
    for i in range(1, N):
        nR[i] = nR[0] + (nR[N] - nR[0]) * i / N
        nE[i] = nE[0] + (nE[N] - nE[0]) * (N - i) / N

    K_profile = np.zeros((N, 3))
    err, n_iter = 1.0, 0

    while err > tol and n_iter < max_iter:
        n_iter += 1
        xR_old, xE_old = xR.copy(), xE.copy()

        for s in range(N):
            gR = nrtl_gamma(xR[s], delta_g, alpha, T, R)
            gE = nrtl_gamma(xE[s + 1], delta_g, alpha, T, R)
            K = gE / gR
            K_profile[s] = K

            # Равновесный рафинат + поправка Мёрфри
            xR_eq = K * xE[s + 1]
            if xR_eq.sum() > 0:
                xR_eq /= xR_eq.sum()
            xR_new = xR[s] + E_stage * (xR_eq - xR[s])
            xR_new = np.maximum(xR_new, 0)
            if xR_new.sum() > 0:
                xR_new /= xR_new.sum()
            xR[s + 1] = damp * xR_new + (1 - damp) * xR[s + 1]
            xR[s + 1] = np.maximum(xR[s + 1], 0)
            xR[s + 1] /= xR[s + 1].sum()

            # Равновесный экстракт
            gR2 = nrtl_gamma(xR[s + 1], delta_g, alpha, T, R)
            gE2 = nrtl_gamma(xE[s], delta_g, alpha, T, R)
            K2 = gR2 / gE2
            xE_eq = xR[s + 1] / K2
            if xE_eq.sum() > 0:
                xE_eq /= xE_eq.sum()
            xE_new = xE[s + 1] + E_stage * (xE_eq - xE[s + 1])
            xE_new = np.maximum(xE_new, 0)
            if xE_new.sum() > 0:
                xE_new /= xE_new.sum()
            xE[s] = damp * xE_new + (1 - damp) * xE[s]
            xE[s] = np.maximum(xE[s], 0)
            xE[s] /= xE[s].sum()

        # Пересчёт потоков из общего баланса
        mat = np.array([[1.0, 1.0], [xE[0, 0], xR[N, 0]]])
        rhs = np.array([n_F + n_S, n_F * x_F[0] + n_S * x_S[0]])
        if abs(np.linalg.det(mat)) > 1e-10:
            sol = np.linalg.solve(mat, rhs)
            nE[0] = max(sol[0], 0.0)
            nR[N] = max(sol[1], 0.0)
        for i in range(N):
            nE[i + 1] = max(nE[0] + nR[0] - nR[i + 1], 0.0)

        err = np.sum(np.abs(xR - xR_old)) + np.sum(np.abs(xE - xE_old))

    # --- Массоперенос и габариты ---
    N_real = int(np.ceil(N / E_stage))
    col_height = N_real * HETS

    C_flood = 0.035
    u_flood = C_flood * (sigma / 30) ** 0.2 * np.sqrt((rho_R - rho_E) / rho_E) \
        if rho_R > rho_E else \
        C_flood * (sigma / 30) ** 0.2 * np.sqrt(abs(rho_R - rho_E) / rho_E)
    u_oper = 0.4 * u_flood

    M_mix_R = float(x_F @ _M)
    Q_vol_R = n_F * M_mix_R / rho_R / 3600.0
    A_col = Q_vol_R / max(u_oper, 1e-6)
    d_col = max(2 * np.sqrt(A_col / np.pi), 0.05)

    D_mol, d_drop = 1.5e-9, 2e-3
    ratio = u_oper / max(u_flood, 1e-6)
    K_La = 0.2 * np.sqrt(D_mol) * ratio ** 0.7 / d_drop ** 2
    NTU = K_La * col_height / max(u_oper, 1e-6)
    HTU = col_height / max(NTU, 1e-6)

    raff = Stream(nR[N], xR[N].copy(), T, "Рафинат")
    extr = Stream(nE[0], xE[0].copy(), T, "Экстракт")

    if verbose:
        print(f"\n--- Экстракционная колонна ({N} теор. ступеней) ---")
        print(f"  Итерации: {n_iter}  |  ошибка: {err:.2e}")
        print(f"  Рафинат:  n={raff.n:.2f} кмоль/ч  "
              f"xA={raff.x[0]:.4f}  xB={raff.x[1]:.5f}  xC={raff.x[2]:.5f}")
        print(f"  Экстракт: n={extr.n:.2f} кмоль/ч  "
              f"xA={extr.x[0]:.4f}  xB={extr.x[1]:.4f}  xC={extr.x[2]:.4f}")
        print(f"  Эффективность Мёрфри: {E_stage*100:.0f}%  =>  реальных ступеней: {N_real}")
        print(f"  Высота: {col_height:.2f} м   Диаметр: {d_col:.3f} м")
        print(f"  Скорость захлёб.: {u_flood:.4f} м/с   Рабочая: {u_oper:.4f} м/с")
        print(f"  K_La={K_La:.2e} 1/с   NTU={NTU:.2f}   HTU={HTU:.4f} м")

    return ExtractorResult(
        raffinate=raff, extract=extr,
        stage_xR=xR, stage_xE=xE, stage_nR=nR, stage_nE=nE,
        K_profile=K_profile, N_stages_req=N_real, E_stage=E_stage,
        col_height=col_height, col_diameter=d_col,
        flooding_vel=u_flood, oper_vel=u_oper,
        K_La=K_La, NTU=NTU, HTU=HTU, HETS=HETS,
        iters=n_iter, error=err,
    )
