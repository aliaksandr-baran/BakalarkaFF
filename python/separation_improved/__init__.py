"""separation_improved — Улучшенная версия пакета моделирования разделения.

Отличия от оригинального separation:
  * extractor_column заменён на extractor_column_newton — одновременный
    Newton/TRF-решатель (scipy.optimize.least_squares) вместо итерации
    с демпфированием. Устраняет численную неустойчивость MATLAB-порта.
  * bubble_T заменён на bubble_T_improved — метод Брента (scipy.optimize.brentq)
    вместо ньютоновского; надёжен при низких давлениях (≈2.45 кПа).
  * vle_flash заменён на vle_flash_improved — использует улучшенный bubble_T.

Все остальные подсистемы (NRTL, LLE, насос, теплообменник, испаритель,
дроссельный клапан, экономика) импортируются из оригинального separation.
"""

# Оригинальные подсистемы (без изменений)
from separation.nrtl import nrtl_gamma, R_GAS
from separation.lle import lle_solver, LLEResult
from separation.extractor import Stream, ExtractorResult
from separation.mass_transfer import mass_transfer_analysis, MassTransferResult
from separation.valve import throttle_valve, ValveResult
from separation.evaporator import evaporator_calc, EvaporatorResult
from separation.heat_exchanger import heat_exchanger, HXResult
from separation.pump import pump_sizing, PumpResult
from separation.economics import process_economics, Equipment, UtilityRates, EconomicsResult
from separation.vle import psat, FlashResult
from separation import thermo

# Улучшенные подсистемы
from .extractor_newton import extractor_column_newton as extractor_column
from .vle_improved import bubble_T_improved as bubble_T, vle_flash_improved as vle_flash

__all__ = [
    # NRTL
    "nrtl_gamma", "R_GAS",
    # LLE
    "lle_solver", "LLEResult",
    # VLE (улучшенные)
    "vle_flash", "bubble_T", "psat", "FlashResult",
    # Экстрактор (Newton)
    "extractor_column", "Stream", "ExtractorResult",
    # Прочие подсистемы (оригинальные)
    "mass_transfer_analysis", "MassTransferResult",
    "throttle_valve", "ValveResult",
    "evaporator_calc", "EvaporatorResult",
    "heat_exchanger", "HXResult",
    "pump_sizing", "PumpResult",
    "process_economics", "Equipment", "UtilityRates", "EconomicsResult",
    "thermo",
]
