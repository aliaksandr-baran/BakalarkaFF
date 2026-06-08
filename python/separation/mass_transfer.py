"""mass_transfer.py — Детальный анализ массопереноса в экстракторе.

Порт MATLAB-функции ``matlab/functions/mass_transfer_analysis.m``.
Рассчитывает K, движущие силы, эффективность Мёрфри, фактор экстракции,
селективность по каждой ступени.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .nrtl import nrtl_gamma, R_GAS


@dataclass
class MassTransferResult:
    K: np.ndarray
    driving_force_R: np.ndarray
    E_Murphree_R: np.ndarray
    lambda_factor: np.ndarray
    NTU_per_stage: np.ndarray
    NTU_total: float
    selectivity: np.ndarray


def mass_transfer_analysis(
    stage_xR: np.ndarray,
    stage_xE: np.ndarray,
    stage_nR: np.ndarray,
    stage_nE: np.ndarray,
    delta_g: np.ndarray,
    alpha: float,
    T: float,
    R: float = R_GAS,
    *,
    HETS: float = 0.15,
    verbose: bool = True,
) -> MassTransferResult:
    """Постадийный анализ массопереноса.

    K_i = gamma_E_i / gamma_R_i
    Движущая сила (рафинат): x - x*
    EMV = (x_in - x_out) / (x_in - x*_out)
    Фактор экстракции: lambda_i = K_i * (L_E/L_R)
    Селективность: beta_BA = K_B / K_A
    """
    N = stage_xR.shape[0] - 1

    K_mat = np.zeros((N, 3))
    DF_R = np.zeros((N, 3))
    E_MV_R = np.zeros((N, 3))
    lam = np.zeros((N, 3))
    NTU_stage = np.zeros(N)
    sel = np.zeros(N)

    if verbose:
        print("\n--- Анализ массопереноса по ступеням экстрактора ---")
        print(f"{'Ст':<3} {'K_A':>7} {'K_B':>8} {'K_C':>9} "
              f"{'ΔxA':>9} {'ΔxB':>9} {'EMV_A':>7} {'EMV_B':>7} {'λ_B':>6} {'β_BA':>6}")
        print("-" * 80)

    for s in range(N):
        xR_in = stage_xR[s]
        xR_out = stage_xR[s + 1]
        xE_in = stage_xE[s + 1]
        xE_out = stage_xE[s]

        gR = nrtl_gamma(xR_out, delta_g, alpha, T, R)
        gE = nrtl_gamma(xE_out, delta_g, alpha, T, R)
        K = gE / gR
        K_mat[s] = K

        xR_eq = K * xE_out
        if xR_eq.sum() > 0:
            xR_eq /= xR_eq.sum()
        xE_eq = xR_out / np.maximum(K, 1e-15)
        if xE_eq.sum() > 0:
            xE_eq /= xE_eq.sum()

        DF_R[s] = xR_out - xR_eq

        denom_R = xR_in - xR_eq
        for ci in range(3):
            E_MV_R[s, ci] = ((xR_in[ci] - xR_out[ci]) / denom_R[ci]
                             if abs(denom_R[ci]) > 1e-8 else 1.0)

        L_R = max(stage_nR[s], 1e-9)
        L_E = max(stage_nE[s + 1], 1e-9)
        lam[s] = K * (L_E / L_R)

        avg_df = max(abs(DF_R[s, 1]), 1e-9)
        NTU_stage[s] = abs(xR_in[1] - xR_out[1]) / avg_df

        sel[s] = K[1] / K[0] if abs(K[0]) > 1e-9 else np.nan

        if verbose:
            print(f"{s+1:<3} {K[0]:>7.3f} {K[1]:>8.3f} {K[2]:>9.4f} "
                  f"{DF_R[s,0]:>9.5f} {DF_R[s,1]:>9.5f} "
                  f"{E_MV_R[s,0]:>7.3f} {E_MV_R[s,1]:>7.3f} "
                  f"{lam[s,1]:>6.3f} {sel[s]:>6.2f}")

    NTU_total = float(NTU_stage.sum())

    if verbose:
        print("-" * 80)
        print(f"  Суммарное NTU (по MeOH): {NTU_total:.3f}")
        valid = np.isfinite(E_MV_R[:, 1])
        print(f"  Средняя эффективность Мёрфри EMV_MeOH: "
              f"{np.mean(E_MV_R[valid,1])*100:.1f}%")
        print(f"  Средняя селективность β_BA: {np.nanmean(sel):.2f}")

    return MassTransferResult(
        K=K_mat, driving_force_R=DF_R, E_Murphree_R=E_MV_R,
        lambda_factor=lam, NTU_per_stage=NTU_stage,
        NTU_total=NTU_total, selectivity=sel,
    )
