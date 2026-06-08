"""extractor_newton.py — Улучшенный решатель противоточного ЖЖ-экстрактора.

Вместо итерации с демпфированием (порт MATLAB) используется одновременный
ньютоновский решатель: все N ступеней формулируются как система 6N нелинейных
уравнений (3 баланса + 3 уравнения равновесия на ступень) и решаются
scipy.optimize.fsolve за один вызов.

Результат: ExtractorResult — совместим с оригинальным пакетом separation.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import fsolve, least_squares

from separation.nrtl import nrtl_gamma, R_GAS
from separation.extractor import Stream, ExtractorResult, _M


# ---------------------------------------------------------------------------
# Вспомогательные функции
# ---------------------------------------------------------------------------

def _pack(r: np.ndarray, e: np.ndarray) -> np.ndarray:
    """r (N×3) и e (N×3) → вектор длиной 6N."""
    return np.concatenate([r.ravel(), e.ravel()])


def _unpack(v: np.ndarray, N: int) -> tuple[np.ndarray, np.ndarray]:
    """Вектор длиной 6N → r (N×3), e (N×3)."""
    r = v[:3 * N].reshape(N, 3)
    e = v[3 * N:].reshape(N, 3)
    return r, e


def _residuals(
    v: np.ndarray,
    N: int,
    r0_feed: np.ndarray,
    e_Np1_solvent: np.ndarray,
    delta_g: np.ndarray,
    alpha: float,
    T: float,
    R: float,
    n_ref: float,
) -> np.ndarray:
    """6N масштабированных невязок для системы N ступеней.

    Ступени пронумерованы s=0..N-1 (s=0 — верх, s=N-1 — низ).
    Питание (рафинат) поступает сверху: r_feed = r_{-1}.
    Растворитель (экстракт) поступает снизу: e_solvent = e_{N}.

    Уравнения на ступень s:
      (0-2) Баланс по компонентам (нормированный на n_ref):
              (r[s] + e[s] - r[s-1] - e[s+1]) / n_ref = 0
      (3-5) Равновесие (K_i * xE_i - xR_i = 0):
              K_i = gamma_E_i / gamma_R_i
    """
    r, e = _unpack(v, N)

    res = np.zeros(6 * N)

    for s in range(N):
        r_in = r0_feed if s == 0 else r[s - 1]
        e_in = e_Np1_solvent if s == N - 1 else e[s + 1]

        # --- Материальный баланс (масштаб: n_ref) ---
        res[6 * s: 6 * s + 3] = (r[s] + e[s] - r_in - e_in) / n_ref

        # --- Равновесие K_i * xE_i = xR_i ---
        nR_s = max(r[s].sum(), 1e-14)
        nE_s = max(e[s].sum(), 1e-14)
        xR = np.maximum(r[s] / nR_s, 1e-15)
        xE = np.maximum(e[s] / nE_s, 1e-15)
        xR /= xR.sum()
        xE /= xE.sum()

        gR = nrtl_gamma(xR, delta_g, alpha, T, R)
        gE = nrtl_gamma(xE, delta_g, alpha, T, R)
        K = gE / gR  # K_i = x_R_i / x_E_i (рафинат/экстракт)

        res[6 * s + 3: 6 * s + 6] = K * xE - xR

    return res


def _initial_guess(
    N: int,
    n_F: float,
    x_F: np.ndarray,
    n_S: float,
    x_S: np.ndarray,
) -> np.ndarray:
    """Линейная интерполяция состава вдоль колонны."""
    r = np.zeros((N, 3))
    e = np.zeros((N, 3))

    # Типичный рафинатный профиль: от x_F до ~[0.97, 0.02, 0.01]
    xR_bot = np.array([0.97, 0.02, 0.01])
    # Типичный экстрактный профиль: от ~[0.05, 0.55, 0.40] до x_S
    xE_top = np.array([0.05, 0.55, 0.40])

    # Оценка потоков из общего баланса
    nR_est = n_F * 0.55
    nE_est = n_S + n_F * 0.45

    for s in range(N):
        t = (s + 0.5) / N  # 0 → верх, 1 → низ
        xR_s = (1 - t) * x_F + t * xR_bot
        xR_s = np.maximum(xR_s, 1e-6)
        xR_s /= xR_s.sum()

        xE_s = (1 - t) * xE_top + t * x_S
        xE_s = np.maximum(xE_s, 1e-6)
        xE_s /= xE_s.sum()

        r[s] = nR_est * xR_s
        e[s] = nE_est * xE_s

    return _pack(r, e)


# ---------------------------------------------------------------------------
# Основная функция
# ---------------------------------------------------------------------------

def extractor_column_newton(
    feed: Stream,
    solvent: Stream,
    N_stages: int,
    delta_g: np.ndarray,
    alpha: float,
    T: float,
    R: float = R_GAS,
    *,
    E_stage: float = 0.85,
    HETS: float = 0.15,
    rho_R: float = 750.0,
    rho_E: float = 1230.0,
    sigma: float = 15.0,
    tol: float = 1e-9,
    max_iter: int = 2000,
    verbose: bool = True,
) -> ExtractorResult:
    """Противоточный ЖЖ-экстрактор с Newton-решателем (scipy.optimize).

    Формулирует систему 6·N нелинейных уравнений (материальный баланс + NRTL-
    равновесие на каждую ступень) и решает одновременно методом Ньютона.
    Это устраняет проблему расходимости итерационного метода оригинального порта.

    Параметры
    ---------
    feed, solvent : Stream
        Питание и растворитель.
    N_stages : int
        Число теоретических ступеней.
    delta_g, alpha, T : float
        Параметры NRTL.
    E_stage : float
        Эффективность Мёрфри.
    """
    N = N_stages
    n_F, x_F = feed.n, feed.x.copy()
    n_S, x_S = solvent.n, solvent.x.copy()

    # Граничные условия (потоки компонентов)
    r0_feed = n_F * x_F
    e_Np1_solv = n_S * x_S

    # Начальное приближение
    v0 = _initial_guess(N, n_F, x_F, n_S, x_S)

    # Нижняя граница: потоки >= 0 (используем least_squares для bounds)
    lb = np.zeros(6 * N)
    ub = np.full(6 * N, (n_F + n_S) * 1.5)

    n_ref = n_F + n_S  # масштабирующий поток

    def residuals_wrapped(v):
        return _residuals(v, N, r0_feed, e_Np1_solv, delta_g, alpha, T, R, n_ref)

    result = least_squares(
        residuals_wrapped,
        v0,
        bounds=(lb, ub),
        method="trf",
        ftol=tol,
        xtol=tol,
        gtol=tol,
        max_nfev=max_iter * 6 * N,
        verbose=0,
    )

    r_sol, e_sol = _unpack(result.x, N)

    # Потоки на выходе из колонны
    # Рафинат выходит снизу (ступень N-1 → bottom raffinate = r[N-1])
    # Экстракт выходит сверху (ступень 0 → top extract = e[0])
    r_out = np.maximum(r_sol[N - 1], 0.0)
    e_out = np.maximum(e_sol[0], 0.0)

    nR_out = r_out.sum()
    nE_out = e_out.sum()

    xR_out = r_out / max(nR_out, 1e-14)
    xE_out = e_out / max(nE_out, 1e-14)

    # Профили для вывода (N+1 строк, как в оригинале)
    stage_xR = np.vstack([x_F, r_sol / np.maximum(r_sol.sum(axis=1, keepdims=True), 1e-14)])
    stage_xE = np.vstack([e_sol / np.maximum(e_sol.sum(axis=1, keepdims=True), 1e-14), x_S])
    stage_nR = np.concatenate([[n_F], r_sol.sum(axis=1)])
    stage_nE = np.concatenate([e_sol.sum(axis=1), [n_S]])

    # Коэффициенты распределения K = xR/xE по ступеням
    K_profile = np.zeros((N, 3))
    for s in range(N):
        xR_s = stage_xR[s + 1]
        xE_s = stage_xE[s]
        K_profile[s] = np.where(xE_s > 1e-10, xR_s / xE_s, 0.0)

    # --- Массоперенос и габариты ---
    N_real = int(np.ceil(N / E_stage))
    col_height = N_real * HETS

    C_flood = 0.035
    drho = abs(rho_R - rho_E)
    u_flood = C_flood * (sigma / 30) ** 0.2 * np.sqrt(drho / max(rho_E, 1.0))
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

    norm_res = float(np.linalg.norm(result.fun))

    raff = Stream(nR_out, xR_out, T, "Рафинат")
    extr = Stream(nE_out, xE_out, T, "Экстракт")

    if verbose:
        print(f"\n--- Экстракционная колонна Newton ({N} теор. ступеней) ---")
        print(f"  Статус решателя: {result.message}")
        print(f"  Норма невязки: {norm_res:.2e}  |  nfev={result.nfev}")
        print(f"  Рафинат:  n={raff.n:.4f} кмоль/ч  "
              f"xA={raff.x[0]:.4f}  xB={raff.x[1]:.5f}  xC={raff.x[2]:.5f}")
        print(f"  Экстракт: n={extr.n:.4f} кмоль/ч  "
              f"xA={extr.x[0]:.4f}  xB={extr.x[1]:.4f}  xC={extr.x[2]:.4f}")
        print(f"  Эффективность Мёрфри: {E_stage*100:.0f}%  =>  реальных ступеней: {N_real}")
        print(f"  Высота: {col_height:.2f} м   Диаметр: {d_col:.3f} м")
        print(f"  Скорость захлёб.: {u_flood:.4f} м/с   Рабочая: {u_oper:.4f} м/с")
        print(f"  K_La={K_La:.2e} 1/с   NTU={NTU:.2f}   HTU={HTU:.4f} м")

    return ExtractorResult(
        raffinate=raff, extract=extr,
        stage_xR=stage_xR, stage_xE=stage_xE,
        stage_nR=stage_nR, stage_nE=stage_nE,
        K_profile=K_profile, N_stages_req=N_real, E_stage=E_stage,
        col_height=col_height, col_diameter=d_col,
        flooding_vel=u_flood, oper_vel=u_oper,
        K_La=K_La, NTU=NTU, HTU=HTU, HETS=HETS,
        iters=result.nfev, error=norm_res,
    )
