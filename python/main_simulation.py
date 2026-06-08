#!/usr/bin/env python3
"""main_simulation.py — Главный скрипт моделирования.

Разделение MTBE–Метанол с помощью ионной жидкости [BMIM][HSO4].
Порт MATLAB-скрипта ``matlab/main_simulation.m``.

Запуск:
    python main_simulation.py

Структура:
    1. Параметры системы (импорт из separation.thermo)
    2. Экстрактор (ЖЖЭ)
    3. Анализ массопереноса
    4. Клапан V-1
    5. ЖПЭ-расчёт паров
    6. Испаритель D-1 (MeOH + рег. ИЖ)
    7. Испаритель D-2 (MTBE)
    8. Теплообменники
    9. Насосы
   10. Экономика
   11. Сводная таблица потоков
   12. Итоговые показатели
"""
from __future__ import annotations

import numpy as np

import separation as sep
from separation import thermo


def main() -> None:
    print("=" * 57)
    print(" МОДЕЛИРОВАНИЕ: Разделение MTBE-MeOH с [BMIM][HSO4]")
    print(" Python-версия с расчётом массо- и теплопереноса")
    print("=" * 57)

    # --- РАЗДЕЛ 1: Параметры ---
    R = sep.R_GAS
    T_ext = 318.15
    alpha = thermo.ALPHA
    delta_g = thermo.DELTA_G
    ant = thermo.ANT

    n_F, x_F = 531.1327, np.array([0.54, 0.46, 0.00])
    n_S, x_S = 258.3784, np.array([0.001559, 0.019959, 0.978483])
    hrs_year = 7920

    # --- РАЗДЕЛ 2: Экстрактор ---
    print("\n=== 2. Экстрактор (5 ступеней, ЖЖЭ) ===")
    feed = sep.Stream(n_F, x_F, T_ext, "Питание F")
    solvent = sep.Stream(n_S, x_S, T_ext, "Растворитель S")
    ext = sep.extractor_column(
        feed, solvent, 5, delta_g, alpha, T_ext, R,
        E_stage=0.85, HETS=0.15, rho_R=740, rho_E=1230, sigma=15,
    )

    # Верифицированные данные из Excel
    n6, x6 = 619.130, np.array([0.03619, 0.55568, 0.40813])
    n18, x18 = 287.206, np.array([0.99679, 0.00047, 0.00275])

    # --- РАЗДЕЛ 3: Массоперенос ---
    print("\n=== 3. Анализ массопереноса ===")
    ext_profile = np.array([
        [29.142, 0.03619, 0.55568, 0.40813, 15.778, 0.86352, 0.13521, 0.00127],
        [14.396, 0.01037, 0.16247, 0.82716, 13.854, 0.97561, 0.02383, 0.00056],
        [12.472, 0.00327, 0.04295, 0.95378, 13.593, 0.99367, 0.00585, 0.00048],
        [12.211, 0.00259, 0.02335, 0.97406, 13.554, 0.99638, 0.00315, 0.00047],
        [12.172, 0.00249, 0.02040, 0.97711, 13.519, 0.99679, 0.00275, 0.00047],
    ])
    stage_xE = np.vstack([ext_profile[:, 1:4], x_S])
    stage_xR = np.vstack([x_F, ext_profile[:, 5:8]])
    stage_nE = np.append(ext_profile[:, 0], n_S)
    stage_nR = np.append(n_F, ext_profile[:, 4])

    mt = sep.mass_transfer_analysis(
        stage_xR, stage_xE, stage_nR, stage_nE, delta_g, alpha, T_ext, R,
    )

    # --- РАЗДЕЛ 4: Клапан V-1 ---
    print("\n=== 4. Клапан V-1 (101.3 → 50 кПа) ===")
    sep.throttle_valve(
        n6, x6, 68 + 273.15, 101325, 50e3, ant, delta_g, alpha, R, name="V-1",
    )
    n9, x9 = 116.825, np.array([0.18384, 0.81616, 0.0])
    n11, x11 = n9, x9

    # --- РАЗДЕЛ 5: ЖПЭ паров ---
    print("\n=== 5. ЖПЭ-расчёт паров (поток 9) ===")
    fr9 = sep.vle_flash(x9, 50e3, 45.114 + 273.15, ant, delta_g, alpha, R)
    print(f"  Поток 9: V_frac={fr9.V_frac:.4f}  "
          f"T_пузырьк={fr9.T_bub-273.15:.2f}°C  сошлось={fr9.converged}")
    print(f"  Состав пара: yA={fr9.y[0]:.5f}  yB={fr9.y[1]:.5f}")

    # --- РАЗДЕЛ 6: Испаритель D-1 ---
    print("\n=== 6. Испаритель D-1 (P=2.45 кПа) ===")
    n12, x12 = 502.296, np.array([0.00185, 0.50306, 0.49509])
    sep.evaporator_calc(
        n12, x12, 45.114 + 273.15, 2.45e3, ant, delta_g, alpha, R,
        T_steam=145 + 273.15, U_evap=600, name="D-1",
    )
    n13, x13 = 244.455, np.array([0.003810, 0.996190, 0.0])
    n16, x16 = 257.842, np.array([2.34e-7, 0.020000, 0.980000])

    # --- РАЗДЕЛ 7: Испаритель D-2 ---
    print("\n=== 7. Испаритель D-2 (P=73.5 кПа) ===")
    sep.evaporator_calc(
        n18, x18, 45 + 273.15, 73.5e3, ant, delta_g, alpha, R,
        T_steam=145 + 273.15, U_evap=900, name="D-2",
    )
    n20, x20 = 286.665, np.array([0.99725, 0.002750, 0.0])
    n23w, x23w = 0.537, np.array([0.7497, 0.000314, 0.250])

    # --- РАЗДЕЛ 8: Теплообменники ---
    print("\n=== 8. Теплообменники (LMTD-метод) ===")
    sep.heat_exchanger(145, 145, 25, 45, 437.02, U=750, name="Подогрев H-1")
    sep.heat_exchanger(145, 145, 45, 68, 715.87, U=890, name="Подогрев H-2")
    sep.heat_exchanger(45.114, 37.646, 20, 35, 1159.62, U=1500, name="Конденс K-1")
    sep.heat_exchanger(45.124, 25, 5, 20, 2584.12, U=2250, name="Конденс K-2")
    sep.heat_exchanger(45.044, 25, 5, 20, 2576.54, U=750, name="Конденс K-3")
    sep.heat_exchanger(145, 145, 20, 45.124, 2469.61, U=450, name="Испарит D-1")
    sep.heat_exchanger(145, 145, 20, 45.0, 2276.68, U=900, name="Испарит D-2")

    # --- РАЗДЕЛ 9: Насосы ---
    print("\n=== 9. Расчёт насосов ===")
    sep.pump_sizing(n_F, x_F, 101325, 101325 + 740*9.81*5, dz=5, name="P-1 (питание)")
    sep.pump_sizing(n_S, x_S, 101325, 101325 + 1275*9.81*1, dz=1, name="P-2 (раств.)")
    sep.pump_sizing(n18, x18, 101325, 73500 + 740*9.81*1, dz=1, name="P-3 (рафинат)")
    sep.pump_sizing(n6, x6, 101325, 101325 + 1100*9.81*5, dz=5, name="P-4 (экстракт)")
    sep.pump_sizing(n16, x16, 2450, 101325 + 1275*9.81*0.5, dz=0.5, name="P-5 (рег.ИЖ)")
    sep.pump_sizing(n23w, x23w, 73500, 101325 + 1100*9.81*2, dz=2, name="P-6 (ост.D-2)")

    # --- РАЗДЕЛ 10: Экономика ---
    print("\n=== 10. Экономический расчёт ===")
    names = ['Насос P-1', 'Насос P-2', 'Насос P-3', 'Насос P-4', 'Насос P-5',
             'Насос P-6', 'Экстрактор', 'Подогрев H-1', 'Подогрев H-2',
             'Конденс K-1', 'Конденс K-2', 'Конденс K-3', 'Испарит D-1', 'Испарит D-2']
    base = [3000, 2100, 2170, 9720, 2737, 1771, 43243, 2124, 2479,
            9912, 23270, 70855, 169999, 119029]
    eq_list = [sep.Equipment(name=n, cost_base=c) for n, c in zip(names, base)]
    util = sep.UtilityRates(
        steam_price=12, cooling_price=0.693, elec_price=0.172,
        steam_kg_h=2.80, cooling_kg_h=73.67, elec_kW=4.96,
    )
    eco = sep.process_economics(eq_list, util, hrs_year, IL_cost=1164891)

    # --- РАЗДЕЛ 11: Сводная таблица потоков ---
    n15, x15 = n13, x13
    n22, x22 = n20, x20
    n23m = np.array([0.47601, 0.52378, 0.00021])
    n23n = n_F + n11

    print("\n" + "=" * 58)
    print(" ТАБЛИЦА: СВОДНАЯ ТАБЛИЦА ВСЕХ ПОТОКОВ ПРОЦЕССА")
    print("=" * 58)
    print(f"{'Поток':<7}{'Фаза':<5}{'T[°C]':>8}{'P[кПа]':>8}"
          f"{'n[кмол/ч]':>11}{'xA':>9}{'xB':>9}{'xC':>9}")
    print("-" * 66)
    rows = [
        ("1", "L", 25.0, 101.3, n_F, x_F),
        ("2,3", "L", 45.0, 101.3, n23n, n23m),
        ("5", "L", 45.0, 101.3, n_S, x_S),
        ("6", "L", 45.0, 101.3, n6, x6),
        ("7", "L", 68.0, 101.3, n6, x6),
        ("9", "G", 45.11, 50.0, n9, x9),
        ("11", "L", 37.65, 101.3, n11, x11),
        ("12", "L", 45.11, 50.0, n12, x12),
        ("13", "G", 45.12, 2.45, n13, x13),
        ("15", "L", 25.0, 101.3, n15, x15),
        ("16", "L", 45.12, 2.45, n16, x16),
        ("18", "L", 45.0, 101.3, n18, x18),
        ("20", "G", 45.04, 73.5, n20, x20),
        ("22", "L", 25.0, 101.3, n22, x22),
        ("23", "L", 45.04, 73.5, n23w, x23w),
    ]
    for sid, ph, T, P, n, x in rows:
        print(f"{sid:<7}{ph:<5}{T:>8.2f}{P:>8.2f}{n:>11.3f}"
              f"{x[0]:>9.5f}{x[1]:>9.5f}{x[2]:>9.5f}")

    # --- РАЗДЕЛ 12: Итоговые показатели ---
    recovery_B = (n15 * x15[1]) / (n_F * x_F[1]) * 100
    recovery_A = (n22 * x22[0]) / (n_F * x_F[0]) * 100
    beta_BA = (x6[1] / max(x6[0], 1e-9)) / (x18[1] / max(x18[0], 1e-9))

    print("\n" + "=" * 58)
    print(" ИТОГОВЫЕ ПОКАЗАТЕЛИ ПРОЦЕССА")
    print("=" * 58)
    print(f"  Число теор. ступеней экстрактора:       5")
    print(f"  Число реальных ступеней (EMV=85%):      {ext.N_stages_req}")
    print(f"  Высота экстракц. колонны:               {ext.col_height:.2f} м")
    print(f"  Диаметр экстракц. колонны:              {ext.col_diameter:.3f} м")
    print(f"  Соотношение S/F:                        {n_S/n_F:.4f}")
    print(f"  Степень извлечения MeOH:                {recovery_B:.2f} %")
    print(f"  Степень извлечения MTBE:                {recovery_A:.2f} %")
    print(f"  Селективность β_BA:                     {beta_BA:.2f}")
    print(f"  NTU (суммарное, MeOH):                  {mt.NTU_total:.3f}")
    print(f"  Чистота MeOH (пот.15) xB:               {x15[1]:.5f}")
    print(f"  Чистота MTBE (пот.22) xA:               {x22[0]:.5f}")
    print(f"  Чистота рег. ИЖ (пот.16) xC:           {x16[2]:.5f}")
    print("-" * 57)
    print(f"  CAPEX оборудование:                     {eco.CAPEX_equip:.0f} €")
    print(f"  TCI (полные капиталовложения):          {eco.TCI:.0f} €")
    print(f"  OPEX:                                   {eco.OPEX:.0f} €/год")
    print("=" * 58 + "\n")


if __name__ == "__main__":
    main()
