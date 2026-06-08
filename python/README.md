# Python port

Python re-implementation of the separation model — the final step of the
**Excel → MATLAB → Python** migration.

Faithful port of the modular MATLAB library in [`../matlab/functions`](../matlab/functions),
with all comments in Russian.

## Structure

```
python/
├── separation/              Пакет моделирования
│   ├── __init__.py
│   ├── nrtl.py              Активностные коэффициенты NRTL
│   ├── lle.py               Жидкость-жидкость равновесие (ЖЖЭ)
│   ├── vle.py               Пар-жидкость равновесие (flash, Антуан)
│   ├── extractor.py         Противоточный экстрактор + габариты колонны
│   ├── mass_transfer.py     Анализ массопереноса (K, EMV, λ, селективность)
│   ├── valve.py             Дросселирующий клапан (изоэнтальпийный, Джоуль-Томсон)
│   ├── evaporator.py        Испаритель (тепловой баланс, площадь)
│   ├── heat_exchanger.py    Теплообменник (LMTD, число трубок)
│   ├── pump.py              Насос (гидравлика, мощность, NPSH)
│   ├── economics.py         CAPEX (CEPCI) / OPEX / TCI
│   └── thermo.py            Физические свойства и параметры NRTL/Антуана
├── tests/
│   └── test_separation.py   Юнит-тесты (pytest)
├── main_simulation.py       Главный скрипт (оркестратор всех подсистем)
└── requirements.txt
```

## Usage

```bash
pip install -r requirements.txt

# Полная симуляция процесса
python main_simulation.py

# Тесты
python -m pytest tests/ -v
```

## Mapping MATLAB → Python

| MATLAB (`matlab/functions/`) | Python (`separation/`)        |
|------------------------------|-------------------------------|
| `nrtl_gamma.m`               | `nrtl.nrtl_gamma`             |
| `lle_solver.m`               | `lle.lle_solver`              |
| `vle_flash.m`                | `vle.vle_flash`, `vle.bubble_T` |
| `extractor_column.m`         | `extractor.extractor_column`  |
| `mass_transfer_analysis.m`   | `mass_transfer.mass_transfer_analysis` |
| `throttle_valve.m`           | `valve.throttle_valve`        |
| `evaporator_calc.m`          | `evaporator.evaporator_calc`  |
| `heat_exchanger.m`           | `heat_exchanger.heat_exchanger` |
| `pump_sizing.m`              | `pump.pump_sizing`            |
| `process_economics.m`        | `economics.process_economics` |
| `main_simulation.m`          | `main_simulation.py`          |

## Model notes

- **NRTL** — interaction parameters `delta_g` [kJ/mol], non-randomness `alpha = 0.4`.
- **LLE** — successive substitution on `K = gamma_E / gamma_R` with damping.
- **VLE** — Rachford-Rice flash + Antoine `Psat`; ionic liquid (C) is non-volatile.
- Stream table uses the Excel-verified compositions; the equilibrium solvers
  reproduce the MATLAB algorithms (the extractor and low-pressure bubble-point
  routines carry the same numerical behaviour as their MATLAB source).
