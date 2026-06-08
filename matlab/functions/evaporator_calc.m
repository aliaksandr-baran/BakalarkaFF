function result = evaporator_calc(feed, P_op, ant, delta_g, alpha, R, opts)
%EVAPORATOR_CALC  Расчёт испарителя (однократное испарение / ЖПЭ)
%                 Совместное решение материального и теплового балансов
%
%  Испаритель реализуется как однократное равновесное испарение при
%  заданном давлении; температура кипения определяется итерационно.
%
%  Входные параметры:
%    feed   — структура входного потока:
%               .n  — молярный поток [кмоль/ч]
%               .x  — мольные доли [xA, xB, xC]
%               .T  — температура подачи [К]
%               .name — имя потока
%    P_op   — рабочее давление испарителя [Па]
%    ant    — матрица Антуана [2×3] (строки: A=MTBE, B=MeOH)
%    delta_g — матрица NRTL [кДж/моль]
%    alpha   — параметр нерандомности NRTL
%    R       — газовая постоянная [Дж/моль/К]
%    opts    — (необязательно):
%               .T_steam    — температура греющего пара [К] (умолч. 418)
%               .V_spec     — заданная доля пара (если NaN — рассчитывается из P)
%               .cp_coeff   — [cpA_kgK, cp_B_dippr_coeff(1×4), cp_IL_coeff(1×3)]
%               .vap_A      — структура параметров теплоты исп. MTBE
%               .vap_B      — структура параметров теплоты исп. MeOH
%               .MA, .MB, .MC — мол. массы [кг/кмоль]
%               .U_evap     — коэфф. теплоотдачи испарителя [Вт/м²/К] (600)
%
%  Выходные параметры:
%    result — структура:
%      .vapor     — поток пара:    .n, .x, .T
%      .liquid    — поток жидкости: .n, .x, .T
%      .T_op      — рабочая температура [К]
%      .V_frac    — мольная доля пара
%      .Q_duty    — тепловая нагрузка [кВт]
%      .Q_sensible — теплота нагрева жидкости [кВт]
%      .Q_latent   — теплота испарения [кВт]
%      .A_evap    — площадь поверхности теплообмена [м²]
%      .LMTD      — LMTD [К]

    if nargin < 7, opts = struct(); end
    T_steam = getf(opts,'T_steam', 418.15);   % 145°C, пар 4 бар
    U_evap  = getf(opts,'U_evap',  600);      % Вт/м²/К (кипение/пар)

    % Молярные массы
    MA = getf(opts,'MA', 88.15);
    MB = getf(opts,'MB', 32.04);
    MC = getf(opts,'MC', 236.29);

    % Параметры теплоёмкости
    cpA_kgK = getf(opts,'cpA_kgK', 2.127);
    cpB_coeff = getf(opts,'cpB_coeff', [0.8382,-0.003231,8.296e-6,-1.689e-10]);
    cpIL_coeff = getf(opts,'cpIL_coeff', [753.28,-3.4195,0.0063]);

    % Параметры теплоты испарения
    vap_A_def.A=46.23; vap_A_def.beta=0.2893; vap_A_def.Tc=497.1;
    vap_B_def.A=45.3; vap_B_def.alpha=-0.31; vap_B_def.beta=0.4241; vap_B_def.Tc=512.6;
    vap_A = getf(opts,'vap_A', vap_A_def);
    vap_B = getf(opts,'vap_B', vap_B_def);

    n_in = feed.n;
    x_in = max(feed.x(:)', 0); x_in = x_in/sum(x_in);
    T_in = feed.T;

    % ----------------------------------------------------------------
    % ШАГ 1: Равновесный расчёт (температура кипения и V_frac)
    % ----------------------------------------------------------------
    V_spec = getf(opts,'V_spec', NaN);
    if isnan(V_spec)
        flash_opts.tol = 1e-9; flash_opts.max_iter = 500;
        [V_frac, y_out, x_out, T_op, ~] = vle_flash(x_in, P_op, ...
            310, ant, delta_g, alpha, R, flash_opts);
        % Уточнение T_op из bubble_point при данном P
        T_op = bubble_T_simple(x_in, P_op, ant, delta_g, alpha, R);
        % Пересчёт V_frac при T_op
        [V_frac, y_out, x_out, ~, ~] = vle_flash(x_in, P_op, T_op, ...
            ant, delta_g, alpha, R, flash_opts);
    else
        V_frac = V_spec;
        T_op = bubble_T_simple(x_in, P_op, ant, delta_g, alpha, R);
        flash_opts.tol = 1e-9;
        [~, y_out, x_out, ~, ~] = vle_flash(x_in, P_op, T_op, ...
            ant, delta_g, alpha, R, flash_opts);
    end

    V_frac = min(max(V_frac, 0), 1);

    % ----------------------------------------------------------------
    % ШАГ 2: Потоки пара и жидкости
    % ----------------------------------------------------------------
    n_vap = n_in * V_frac;
    n_liq = n_in * (1 - V_frac);

    % ----------------------------------------------------------------
    % ШАГ 3: Тепловой баланс
    % ----------------------------------------------------------------
    % Теплоёмкости [кДж/кмоль/К]
    cpA = @(T) cpA_kgK * MA;
    cpB = @(T) (cpB_coeff(1)+cpB_coeff(2)*T+cpB_coeff(3)*T^2+cpB_coeff(4)*T^3)*MB*4.184;
    cpIL= @(T) cpIL_coeff(1)+cpIL_coeff(2)*(T-273.15)+cpIL_coeff(3)*(T-273.15)^2;

    % Среднесмесевая теплоёмкость жидкой фазы [кДж/кмоль/К] при T_in
    cp_mix_in = x_in(1)*cpA(T_in) + x_in(2)*cpB(T_in) + x_in(3)*cpIL(T_in);

    % Теплота нагрева жидкости от T_in до T_op [кДж/кмоль]
    if abs(T_op - T_in) > 0.1
        dH_sensible = integral(@(T) x_in(1)*cpA(T)+x_in(2)*cpB(T)+x_in(3)*cpIL(T), ...
            T_in, T_op);
    else
        dH_sensible = 0;
    end
    Q_sensible = n_in * dH_sensible / 3600;   % кВт (кДж/ч → кВт)

    % Теплоты испарения при T_op [кДж/моль]
    Hvap_A = vap_A.A * max((vap_A.Tc-T_op)/(vap_A.Tc-298.15), 0)^vap_A.beta;
    exponent_B = vap_B.alpha + vap_B.beta*(T_op/vap_B.Tc);
    Hvap_B = vap_B.A * max((vap_B.Tc-T_op)/(vap_B.Tc-298.15), 0)^exponent_B;

    % Теплота испарения смеси (взвешенная по пару) [кДж/кмоль]
    if sum(y_out(1:2)) > 1e-9
        Hvap_mix = (y_out(1)*Hvap_A + y_out(2)*Hvap_B) * 1000;   % → кДж/кмоль
    else
        Hvap_mix = 0;
    end

    Q_latent = n_vap * Hvap_mix / 3600;   % кВт

    Q_total = Q_sensible + Q_latent;

    % ----------------------------------------------------------------
    % ШАГ 4: Площадь теплообмена
    % ----------------------------------------------------------------
    LMTD_evap = lmtd(T_steam-273.15, T_steam-273.15, T_in-273.15, T_op-273.15);
    if LMTD_evap > 0
        A_evap = Q_total * 1000 / (U_evap * LMTD_evap);
    else
        A_evap = NaN;
    end

    % ----------------------------------------------------------------
    % Результат
    % ----------------------------------------------------------------
    result.vapor.n   = n_vap;
    result.vapor.x   = y_out;
    result.vapor.T   = T_op;
    result.liquid.n  = n_liq;
    result.liquid.x  = x_out;
    result.liquid.T  = T_op;
    result.T_op      = T_op;
    result.V_frac    = V_frac;
    result.Q_duty    = Q_total;
    result.Q_sensible = Q_sensible;
    result.Q_latent  = Q_latent;
    result.A_evap    = A_evap;
    result.LMTD      = LMTD_evap;
    result.Hvap_A    = Hvap_A;
    result.Hvap_B    = Hvap_B;

    % --- Печать ---
    name = getf(feed,'name','Испаритель');
    fprintf('  %-20s  P=%6.2f кПа  T_кип=%6.2f°C  V=%.3f\n', ...
        name, P_op/1e3, T_op-273.15, V_frac);
    fprintf('    Q_сенс=%7.2f кВт  Q_лат=%7.2f кВт  Q_total=%7.2f кВт\n', ...
        Q_sensible, Q_latent, Q_total);
    fprintf('    Площадь: %.2f м²   LMTD: %.2f К\n', A_evap, LMTD_evap);
    fprintf('    Пар:  n=%.2f кмоль/ч  xA=%.4f  xB=%.4f  xC=%.4f\n', ...
        n_vap, y_out(1), y_out(2), y_out(3));
    fprintf('    Жидк: n=%.2f кмоль/ч  xA=%.4f  xB=%.4f  xC=%.4f\n', ...
        n_liq, x_out(1), x_out(2), x_out(3));
end

% -----------------------------------------------------------------------
function T_b = bubble_T_simple(z, P, ant, delta_g, alpha, R)
    T_b = 320;
    for k = 1:300
        Ps_A = exp(ant(1,1)+ant(1,2)/(T_b+ant(1,3)));
        Ps_B = exp(ant(2,1)+ant(2,2)/(T_b+ant(2,3)));
        gm = nrtl_gamma(z, delta_g, alpha, R, T_b);
        f  = gm(1)*Ps_A/P*z(1) + gm(2)*Ps_B/P*z(2) - 1;
        dT = 0.1;
        Ps_A2=exp(ant(1,1)+ant(1,2)/(T_b+dT+ant(1,3)));
        Ps_B2=exp(ant(2,1)+ant(2,2)/(T_b+dT+ant(2,3)));
        gm2 = nrtl_gamma(z,delta_g,alpha,R,T_b+dT);
        f2 = gm2(1)*Ps_A2/P*z(1)+gm2(2)*Ps_B2/P*z(2)-1;
        df = (f2-f)/dT;
        if abs(df)<1e-15, break; end
        T_n = T_b - f/df;
        if abs(T_n-T_b)<1e-6, T_b=T_n; break; end
        T_b = T_n;
    end
end

function L = lmtd(Th_in,Th_out,Tc_in,Tc_out)
    dT1 = Th_in-Tc_out; dT2 = Th_out-Tc_in;
    if abs(dT1-dT2)<1e-4, L=dT1; return; end
    if dT1<=0||dT2<=0, L=max(abs(dT1),abs(dT2)); return; end
    L = (dT1-dT2)/log(dT1/dT2);
end

function val = getf(s,f,d)
    if isfield(s,f), val=s.(f); else, val=d; end
end
