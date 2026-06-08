function result = heat_exchanger(hot, cold, Q_duty, opts)
%HEAT_EXCHANGER  Тепловой расчёт теплообменника (метод LMTD)
%                Определение площади, числа трубок/пластин, проверка U
%
%  Входные параметры:
%    hot  — структура горячего потока:
%              .T_in  — температура на входе [°C]
%              .T_out — температура на выходе [°C]
%              .name  — имя потока (строка)
%    cold — структура холодного потока (те же поля)
%    Q_duty — тепловая нагрузка [кВт]
%    opts — (необязательно) структура:
%              .type  — 'shell_tube' или 'plate' (по умолч. 'shell_tube')
%              .U     — коэффициент теплопередачи [Вт/м²/К] (по умолч. авто)
%              .passes — число ходов по трубной стороне (1 или 2, умолч. 1)
%              .d_tube — наружный диаметр трубок [м] (0.025)
%              .L_tube — длина трубок [м] (3.0)
%              .fouling_R — суммарное сопротивление загрязнения [м²·К/Вт] (1e-4)
%
%  Выходные параметры:
%    result — структура:
%      .Q         — тепловая нагрузка [кВт]
%      .LMTD      — логарифмический температурный напор [К]
%      .F_corr    — поправочный коэффициент F (многоходовые т/о)
%      .U_design  — расчётный U [Вт/м²/К]
%      .A_req     — требуемая площадь теплообмена [м²]
%      .n_tubes   — число трубок (для кожухотрубного)
%      .n_plates  — число пластин (для пластинчатого)
%      .NTU       — число единиц переноса теплоты
%      .effectiveness — эффективность теплообменника
%
%  Пример (испаритель D-1):
%    hot  = struct('T_in',145,'T_out',145,'name','пар 3 бар');
%    cold = struct('T_in',20,'T_out',45.1,'name','поток 12');
%    res  = heat_exchanger(hot, cold, 2469.61)

    if nargin < 4, opts = struct(); end
    type    = getf(opts,'type',    'shell_tube');
    passes  = getf(opts,'passes',  1);
    d_tube  = getf(opts,'d_tube',  0.025);   % м
    L_tube  = getf(opts,'L_tube',  3.0);     % м
    R_foul  = getf(opts,'fouling_R', 1e-4);  % м²·К/Вт

    Th_in  = hot.T_in;   Th_out = hot.T_out;
    Tc_in  = cold.T_in;  Tc_out = cold.T_out;

    % --- Автоматический выбор U [Вт/м²/К] ---
    if isfield(opts,'U')
        U = opts.U;
    else
        U = select_U(type, Th_in, Tc_in);
    end

    % --- LMTD для противотока ---
    dT1 = Th_in  - Tc_out;   % горячий вход — холодный выход
    dT2 = Th_out - Tc_in;    % горячий выход — холодный вход

    if abs(dT1 - dT2) < 1e-4
        LMTD_cf = dT1;
    elseif dT1 <= 0 || dT2 <= 0
        % Невозможная конфигурация
        LMTD_cf = NaN;
        warning('heat_exchanger: некорректные температуры — пересечение.');
    else
        LMTD_cf = (dT1 - dT2) / log(dT1/dT2);
    end

    % --- Поправочный коэффициент F для многоходовых теплообменников ---
    F = lmtd_correction(Th_in, Th_out, Tc_in, Tc_out, passes);
    LMTD = LMTD_cf * F;

    % --- Эффективный U с учётом загрязнения ---
    U_eff = 1 / (1/U + R_foul);

    % --- Требуемая площадь ---
    Q_W = Q_duty * 1000;   % кВт → Вт
    if isnan(LMTD) || LMTD <= 0
        A_req = NaN;
    else
        A_req = Q_W / (U_eff * LMTD);
    end

    % --- Число трубок (кожухотрубный) ---
    A_tube = pi * d_tube * L_tube;
    n_tubes = ceil(A_req / A_tube);

    % --- Число пластин (пластинчатый, h_plate=0.003 м, b=0.5 м, L=1.5 м) ---
    A_plate = 0.5 * 1.5;   % 0.75 м² на пластину
    n_plates = ceil(A_req / A_plate) + 1;

    % --- NTU и эффективность ---
    % Упрощённо: NTU = U*A/(C_min), эффективность из epsilon-NTU
    % Для конденсаторов/испарителей C_min → 0, epsilon → 1
    if abs(Th_in - Th_out) < 0.5   % конденсация/испарение
        NTU = Inf; effectiveness = 1.0;
    else
        dT_max = max(Th_in - Tc_in, 1e-3);
        effectiveness = (Th_in - Th_out) / dT_max;
        NTU = -log(1 - effectiveness);   % для C_r→0
    end

    % --- Формирование результата ---
    result.Q            = Q_duty;
    result.LMTD_cf      = LMTD_cf;
    result.F_corr       = F;
    result.LMTD         = LMTD;
    result.U_input      = U;
    result.U_design     = U_eff;
    result.A_req        = A_req;
    result.n_tubes      = n_tubes;
    result.n_plates     = n_plates;
    result.NTU          = NTU;
    result.effectiveness = effectiveness;
    result.type         = type;

    % --- Печать ---
    if isfield(hot,'name') && isfield(cold,'name')
        name = sprintf('%s / %s', hot.name, cold.name);
    else
        name = 'Теплообменник';
    end
    fprintf('  %-30s  Q=%7.1f кВт  LMTD=%6.2f K  U=%5.0f Вт/м²/К  A=%6.2f м²', ...
        name, Q_duty, LMTD, U_eff, A_req);
    if strcmp(type,'shell_tube')
        fprintf('  n_труб=%d\n', n_tubes);
    else
        fprintf('  n_пластин=%d\n', n_plates);
    end
end

% -----------------------------------------------------------------------
function U = select_U(type, Th_in, Tc_in)
%SELECT_U  Типовые значения U [Вт/м²/К] по справочнику (Kern, 1950)
    if strcmp(type,'plate')
        U = 3000;   % пластинчатый т/о, жидкость-жидкость
        return
    end
    % Кожухотрубный
    is_cond = false; is_evap = false;
    if abs(Th_in - 145) < 30,  is_cond = true; end   % водяной пар ~3 бар
    if abs(Tc_in - 20)  < 50,  is_evap = false; end

    if is_cond
        U = 1200;   % конденсация пара / жидкость
    else
        U = 800;    % жидкость / жидкость
    end
end

% -----------------------------------------------------------------------
function F = lmtd_correction(Th_in, Th_out, Tc_in, Tc_out, passes)
%LMTD_CORRECTION  Поправочный коэффициент F (TEMA стандарт)
%  F ≈ 1 для противотока и конденсации/испарения.
%  Для однохо-двухходового теплообменника используется аналитическая формула.
    if passes == 1
        F = 1.0; return
    end
    R = (Th_in - Th_out) / max(Tc_out - Tc_in, 1e-6);
    P = (Tc_out - Tc_in)  / max(Th_in - Tc_in, 1e-6);
    if abs(R - 1) < 1e-4
        W = (1-P)/(1 - R*P + 1e-15);
        if W <= 0, F = 1.0; return; end
        S_arg = sqrt(2) * P / (1 - P);
        if S_arg <= 0, F = 1.0; return; end
        F = sqrt(2)*P / ((1-P)*log(max((1+W-P+sqrt(2)*P)/(1+W-P-sqrt(2)*P),1e-6)));
    else
        S = sqrt(R^2+1);
        arg1 = (2/P - 1 - R + S) / (2/P - 1 - R - S);
        if arg1 <= 0 || arg1 == 1, F = 1.0; return; end
        F = S * log((1-P*R)/(max(1-P,1e-9))) / ...
            ((R-1) * log(max(arg1, 1e-9)));
    end
    F = min(max(F, 0.5), 1.0);
end

function val = getf(s,f,d)
    if isfield(s,f), val=s.(f); else, val=d; end
end
