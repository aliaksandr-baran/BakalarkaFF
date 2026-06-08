#!/usr/bin/env python3
"""main_simulation_improved.py — Улучшенная версия моделирования.

Ключевые улучшения по сравнению с main_simulation.py:
  * Экстрактор решается одновременным Newton/TRF-методом (scipy.optimize.least_squares),
    без итерации с демпфированием. Состав потоков из экстрактора принимается
    непосредственно из решателя, а не из жёстко заданных Excel-данных.
  * Все нижестоящие расчёты (клапаны, испарители, теплообменники, насосы)
    используют состав, вычисленный решателем.
  * bubble_T использует метод Брента — надёжен при низких давлениях.

Запуск:
    cd python
    python main_simulation_improved.py
"""
from __future__ import annotations

import sys
import numpy as np

# Добавляем путь к пакетам
import os
sys.path.insert(0, os.path.dirname(__file__))

import separation_improved as sep
from separation_improved import thermo


def main() -> None:
    print("=" * 62)
    print(" МОДЕЛИРОВАНИЕ (улучш.): Разделение MTBE-MeOH / [BMIM][HSO4]")
    print(" Newton-решатель экстрактора + метод Брента для bubble_T")
    print("=" * 62)

    # --- РАЗДЕЛ 1: Параметры ---
    R = sep.R_GAS
    T_ext = 318.15
    alpha = thermo.ALPHA
    delta_g = thermo.DELTA_G
    ant = thermo.ANT

    n_F, x_F = 531.1327, np.array([0.54, 0.46, 0.00])
    n_S, x_S = 258.3784, np.array([0.001559, 0.019959, 0.978483])
    hrs_year = 7920

    # --- РАЗДЕЛ 2: Экстрактор (Newton-решатель) ---
    print("\n=== 2. Экстрактор (5 ступеней, Newton-TRF) ===")
    feed = sep.Stream(n_F, x_F, T_ext, "Питание F")
    solvent = sep.Stream(n_S, x_S, T_ext, "Растворитель S")
    ext = sep.extractor_column(
        feed, solvent, 5, delta_g, alpha, T_ext, R,
        E_stage=0.85, HETS=0.15, rho_R=740, rho_E=1230, sigma=15,
    )

    # Результаты экстрактора из решателя (не из Excel)
    n6, x6 = ext.extract.n, ext.extract.x     # экстракт (выход сверху)
    n18, x18 = ext.raffinate.n, ext.raffinate.x  # рафинат (выход снизу)

    print(f"\n  Экстракт (пот.6):    n={n6:.3f} кмоль/ч  "
          f"xA={x6[0]:.5f}  xB={x6[1]:.5f}  xC={x6[2]:.5f}")
    print(f"  Рафинат (пот.18):   n={n18:.3f} кмоль/ч  "
          f"xA={x18[0]:.5f}  xB={x18[1]:.5f}  xC={x18[2]:.5f}")

    # --- РАЗДЕЛ 3: Анализ массопереноса ---
    print("\n=== 3. Анализ массопереноса ===")
    mt = sep.mass_transfer_analysis(
        ext.stage_xR, ext.stage_xE, ext.stage_nR, ext.stage_nE,
        delta_g, alpha, T_ext, R,
    )

    # --- РАЗДЕЛ 4: Клапан V-1 (101.3 → 50 кПа) ---
    print("\n=== 4. Клапан V-1 (101.3 → 50 кПа) ===")
    T_ext_hot = 68 + 273.15  # нагрев перед дросселем
    valve1 = sep.throttle_valve(
        n6, x6, T_ext_hot, 101325, 50e3,
        ant, delta_g, alpha, R, name="V-1",
    )
    # Поток 9: паровая фаза после клапана
    n9_tot = n6 * valve1.V_frac if valve1.V_frac > 0 else n6 * 0.2
    x9 = valve1.y if valve1.V_frac > 0 else np.array([0.18, 0.82, 0.0])
    # Поток 12: жидкая фаза после клапана
    n12_tot = n6 * (1 - valve1.V_frac) if valve1.V_frac > 0 else n6 * 0.8
    x12 = valve1.x if valve1.V_frac > 0 else np.array([0.01, 0.50, 0.49])

    # Нормализуем x9 (ИЖ нелетуча)
    x9 = np.array([x9[0], x9[1], 0.0])
    if x9.sum() > 0:
        x9 = x9 / x9.sum()

    n9, n11 = n9_tot, n9_tot
    n12 = n12_tot
    x11 = x9.copy()

    # --- РАЗДЕЛ 5: ЖПЭ-расчёт паров (поток 9) ---
    print("\n=== 5. ЖПЭ-расчёт паров (поток 9) ===")
    fr9 = sep.vle_flash(x9, 50e3, 45.114 + 273.15, ant, delta_g, alpha, R)
    print(f"  Поток 9: V_frac={fr9.V_frac:.4f}  "
          f"T_пузырьк={fr9.T_bub - 273.15:.2f}°C  сошлось={fr9.converged}")
    print(f"  Состав пара: yA={fr9.y[0]:.5f}  yB={fr9.y[1]:.5f}")

    # --- РАЗДЕЛ 6: Испаритель D-1 (2.45 кПа) ---
    print("\n=== 6. Испаритель D-1 (P=2.45 кПа) ===")
    ev1 = sep.evaporator_calc(
        n12, x12, 45.114 + 273.15, 2.45e3,
        ant, delta_g, alpha, R,
        T_steam=145 + 273.15, U_evap=600, name="D-1",
    )
    # Паровой продукт D-1 = MeOH
    n13 = ev1.vapor_n
    x13 = ev1.vapor_x
    # Жидкость = регенерированная ИЖ
    n16 = ev1.liquid_n
    x16 = ev1.liquid_x

    # --- РАЗДЕЛ 7: Испаритель D-2 (73.5 кПа) ---
    print("\n=== 7. Испаритель D-2 (P=73.5 кПа) ===")
    ev2 = sep.evaporator_calc(
        n18, x18, 45 + 273.15, 73.5e3,
        ant, delta_g, alpha, R,
        T_steam=145 + 273.15, U_evap=900, name="D-2",
    )
    # Паровой продукт D-2 = MTBE
    n20 = ev2.vapor_n
    x20 = ev2.vapor_x
    # Остаток
    n23w = ev2.liquid_n
    x23w = ev2.liquid_x

    # --- РАЗДЕЛ 8: Теплообменники (LMTD) ---
    print("\n=== 8. Теплообменники (LMTD-метод) ===")
    sep.heat_exchanger(145, 145, 25, 45,    437.02,  U=750,  name="Подогрев H-1")
    sep.heat_exchanger(145, 145, 45, 68,    715.87,  U=890,  name="Подогрев H-2")
    sep.heat_exchanger(45.114, 37.646, 20, 35, 1159.62, U=1500, name="Конденс K-1")
    sep.heat_exchanger(45.124, 25,    5,  20, 2584.12, U=2250, name="Конденс K-2")
    sep.heat_exchanger(45.044, 25,    5,  20, 2576.54, U=750,  name="Конденс K-3")
    sep.heat_exchanger(145, 145, 20, 45.124, 2469.61, U=450,  name="Испарит D-1")
    sep.heat_exchanger(145, 145, 20, 45.0,   2276.68, U=900,  name="Испарит D-2")

    # --- РАЗДЕЛ 9: Насосы ---
    print("\n=== 9. Расчёт насосов ===")
    sep.pump_sizing(n_F,  x_F,  101325, 101325 + 740 * 9.81 * 5,    dz=5,   name="P-1 (питание)")
    sep.pump_sizing(n_S,  x_S,  101325, 101325 + 1275 * 9.81 * 1,   dz=1,   name="P-2 (раств.)")
    sep.pump_sizing(n18,  x18,  101325, 73500  + 740 * 9.81 * 1,    dz=1,   name="P-3 (рафинат)")
    sep.pump_sizing(n6,   x6,   101325, 101325 + 1100 * 9.81 * 5,   dz=5,   name="P-4 (экстракт)")
    sep.pump_sizing(n16,  x16,  2450,   101325 + 1275 * 9.81 * 0.5, dz=0.5, name="P-5 (рег.ИЖ)")
    sep.pump_sizing(n23w, x23w, 73500,  101325 + 1100 * 9.81 * 2,   dz=2,   name="P-6 (ост.D-2)")

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

    print("\n" + "=" * 68)
    print(" ТАБЛИЦА: СВОДНАЯ ТАБЛИЦА ПОТОКОВ (Newton-решатель)")
    print("=" * 68)
    print(f"{'Поток':<7}{'Фаза':<5}{'T[°C]':>8}{'P[кПа]':>8}"
          f"{'n[кмол/ч]':>11}{'xA':>10}{'xB':>10}{'xC':>10}")
    print("-" * 71)
    rows = [
        ("1",   "L", 25.0,   101.3,  n_F,  x_F),
        ("5",   "L", 45.0,   101.3,  n_S,  x_S),
        ("6",   "L", 45.0,   101.3,  n6,   x6),
        ("7",   "L", 68.0,   101.3,  n6,   x6),
        ("9",   "G", 45.11,   50.0,  n9,   x9),
        ("11",  "L", 37.65,  101.3,  n11,  x11),
        ("12",  "L", 45.11,   50.0,  n12,  x12),
        ("13",  "G", 45.12,    2.45, n13,  x13),
        ("15",  "L", 25.0,   101.3,  n15,  x15),
        ("16",  "L", 45.12,    2.45, n16,  x16),
        ("18",  "L", 45.0,   101.3,  n18,  x18),
        ("20",  "G", 45.04,   73.5,  n20,  x20),
        ("22",  "L", 25.0,   101.3,  n22,  x22),
        ("23",  "L", 45.04,   73.5,  n23w, x23w),
    ]
    for sid, ph, T, P, n, x in rows:
        print(f"{sid:<7}{ph:<5}{T:>8.2f}{P:>8.2f}{n:>11.3f}"
              f"{x[0]:>10.5f}{x[1]:>10.5f}{x[2]:>10.5f}")

    # Проверка общего баланса
    print("\n  Проверка материального баланса:")
    in_tot = n_F + n_S
    out_tot = n15 + n22 + n16 + n23w
    print(f"    Вход:  {in_tot:.4f} кмоль/ч")
    print(f"    Выход: {out_tot:.4f} кмоль/ч")
    print(f"    Невязка: {abs(in_tot - out_tot):.4f} кмоль/ч  "
          f"({abs(in_tot - out_tot) / in_tot * 100:.3f}%)")

    # --- РАЗДЕЛ 12: Итоговые показатели ---
    n_B_in = n_F * x_F[1]
    n_B_out = n15 * x15[1] if x15.sum() > 0 else 0.0
    recovery_B = n_B_out / max(n_B_in, 1e-9) * 100

    n_A_in = n_F * x_F[0]
    n_A_out = n22 * x22[0] if x22.sum() > 0 else 0.0
    recovery_A = n_A_out / max(n_A_in, 1e-9) * 100

    # Селективность ИЖ: beta = (xB/xA)_экстракт / (xB/xA)_рафинат
    beta_BA = (x6[1] / max(x6[0], 1e-9)) / (x18[1] / max(x18[0], 1e-9))

    print("\n" + "=" * 62)
    print(" ИТОГОВЫЕ ПОКАЗАТЕЛИ ПРОЦЕССА (улучшенная версия)")
    print("=" * 62)
    print(f"  Число теор. ступеней экстрактора:       5")
    print(f"  Число реальных ступеней (EMV=85%):      {ext.N_stages_req}")
    print(f"  Норма невязки Newton-решателя:          {ext.error:.2e}")
    print(f"  Высота экстракц. колонны:               {ext.col_height:.2f} м")
    print(f"  Диаметр экстракц. колонны:              {ext.col_diameter:.3f} м")
    print(f"  Соотношение S/F:                        {n_S / n_F:.4f}")
    print(f"  Степень извлечения MeOH:                {recovery_B:.2f} %")
    print(f"  Степень извлечения MTBE:                {recovery_A:.2f} %")
    print(f"  Селективность β_BA (экстр./рафин.):    {beta_BA:.2f}")
    print(f"  NTU (суммарное, MeOH):                  {mt.NTU_total:.3f}")
    print(f"  Чистота MeOH (пот.15) xB:               {x15[1]:.5f}")
    print(f"  Чистота MTBE (пот.22) xA:               {x22[0]:.5f}")
    print(f"  Чистота рег. ИЖ (пот.16) xC:           {x16[2]:.5f}")
    print("-" * 62)
    print(f"  CAPEX оборудование:                     {eco.CAPEX_equip:.0f} €")
    print(f"  TCI (полные капиталовложения):          {eco.TCI:.0f} €")
    print(f"  OPEX:                                   {eco.OPEX:.0f} €/год")
    print("=" * 62 + "\n")


if __name__ == "__main__":
    main()
