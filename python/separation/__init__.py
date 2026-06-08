"""separation — Пакет моделирования разделения MTBE–Метанол / [BMIM][HSO4].

Порт MATLAB-библиотеки (matlab/functions) на Python.

Подсистемы:
    nrtl            — активностные коэффициенты NRTL
    lle             — жидкость-жидкость равновесие
    vle             — пар-жидкость равновесие (flash)
    extractor       — противоточный экстрактор
    mass_transfer   — анализ массопереноса
    valve           — дросселирующий клапан (изоэнтальпийный)
    evaporator      — испаритель (тепловой баланс)
    heat_exchanger  — теплообменник (LMTD)
    pump            — насос (гидравлика)
    economics       — CAPEX / OPEX
    thermo          — физические свойства и параметры
"""
from .nrtl import nrtl_gamma, R_GAS
from .lle import lle_solver, LLEResult
from .vle import vle_flash, bubble_T, psat, FlashResult
from .extractor import extractor_column, Stream, ExtractorResult
from .mass_transfer import mass_transfer_analysis, MassTransferResult
from .valve import throttle_valve, ValveResult
from .evaporator import evaporator_calc, EvaporatorResult
from .heat_exchanger import heat_exchanger, HXResult
from .pump import pump_sizing, PumpResult
from .economics import process_economics, Equipment, UtilityRates, EconomicsResult
from . import thermo

__all__ = [
    "nrtl_gamma", "R_GAS",
    "lle_solver", "LLEResult",
    "vle_flash", "bubble_T", "psat", "FlashResult",
    "extractor_column", "Stream", "ExtractorResult",
    "mass_transfer_analysis", "MassTransferResult",
    "throttle_valve", "ValveResult",
    "evaporator_calc", "EvaporatorResult",
    "heat_exchanger", "HXResult",
    "pump_sizing", "PumpResult",
    "process_economics", "Equipment", "UtilityRates", "EconomicsResult",
    "thermo",
]
