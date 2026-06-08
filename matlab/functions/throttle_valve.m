function result = throttle_valve(feed, P_out, ant, delta_g, alpha, R, opts)
%THROTTLE_VALVE  Изоэнтальпийное дросселирование (клапан)
%                Расчёт температуры и фазового состояния после клапана
%
%  Изоэнтальпийный процесс: H_in = H_out
%  Температура после клапана T_out находится из условия равенства энтальпий.
%  При T_out < T_пузырьк(P_out) — двухфазная смесь.
%
%  Входные параметры:
%    feed   — структура:
%               .n  — молярный поток [кмоль/ч]
%               .x  — мольные доли [xA, xB, xC]
%               .T  — температура [К]
%               .P  — давление [Па]
%               .name — имя потока
%    P_out  — давление после клапана [Па]
%    ant    — матрица Антуана [2×3]
%    delta_g — матрица NRTL [кДж/моль]
%    alpha   — параметр NRTL
%    R       — газовая постоянная [Дж/моль/К]
%    opts    — (необязательно):
%               .cpA_kgK    — теплоёмкость MTBE [кДж/кг/К] (2.127)
%               .cpB_coeff  — коэффициенты cp MeOH DIPPR [1×4]
%               .cpIL_coeff — коэффициенты cp ИЖ [1×3]
%               .vap_A, .vap_B — параметры теплоты испарения
%               .MA, .MB, .MC  — молярные массы [кг/кмоль]
%
%  Выходные параметры:
%    result — структура:
%      .T_out   — температура после клапана [К]
%      .P_out   — давление после клапана [Па]
%      .phase   — 'liquid', 'two-phase', или 'vapor'
%      .V_frac  — мольная доля пара
%      .vapor   — поток пара:    .n, .x, .T, .P
%      .liquid  — поток жидкости: .n, .x, .T, .P
%      .dT      — перепад температур [К]
%      .Joule_Thomson — коэффициент Джоуля-Томсона [К/Па]

    if nargin < 7, opts = struct(); end
    cpA_kgK   = getf(opts,'cpA_kgK',   2.127);
    cpB_coeff = getf(opts,'cpB_coeff',  [0.8382,-0.003231,8.296e-6,-1.689e-10]);
    cpIL_coeff= getf(opts,'cpIL_coeff', [753.28,-3.4195,0.0063]);
    MA = getf(opts,'MA',88.15); MB=getf(opts,'MB',32.04); MC=getf(opts,'MC',236.29);

    vap_A_def.A=46.23; vap_A_def.beta=0.2893; vap_A_def.Tc=497.1;
    vap_B_def.A=45.3; vap_B_def.alpha=-0.31; vap_B_def.beta=0.4241; vap_B_def.Tc=512.6;
    vap_A = getf(opts,'vap_A', vap_A_def);
    vap_B = getf(opts,'vap_B', vap_B_def);

    x   = max(feed.x(:)', 0); x = x/sum(x);
    T_in= feed.T;
    P_in= feed.P;
    n   = feed.n;

    cpA  = @(T) cpA_kgK * MA;
    cpB  = @(T) (cpB_coeff(1)+cpB_coeff(2)*T+cpB_coeff(3)*T^2+cpB_coeff(4)*T^3)*MB*4.184;
    cpIL = @(T) cpIL_coeff(1)+cpIL_coeff(2)*(T-273.15)+cpIL_coeff(3)*(T-273.15)^2;
    cp_mix = @(T) x(1)*cpA(T) + x(2)*cpB(T) + x(3)*cpIL(T);

    Hvap_A = @(T) vap_A.A*max((vap_A.Tc-T)/(vap_A.Tc-298.15),0)^vap_A.beta*1000;
    exponent_B = @(T) vap_B.alpha + vap_B.beta*(T/vap_B.Tc);
    Hvap_B = @(T) vap_B.A*max((vap_B.Tc-T)/(vap_B.Tc-298.15),0)^exponent_B(T)*1000;
    Hvap_mix = @(T,y) (y(1)*Hvap_A(T)+y(2)*Hvap_B(T));

    % ------------------------------------------------------------------
    % ШАГ 1: Температура пузырькования при P_out
    % ------------------------------------------------------------------
    T_bub = bubble_T_simple(x, P_out, ant, delta_g, alpha, R);

    % ------------------------------------------------------------------
    % ШАГ 2: Изоэнтальпийный баланс
    % Предположение: жидкость, T_out = T_in - ΔP/cp*Joule-Thomson
    % Точный расчёт: H(T_out, P_out, V_frac) = H(T_in, P_in, 0)
    %
    % Если T_bub < T_in → при P_out уже часть испарится
    % Задача: найти T_out и V_frac из системы:
    %   cp_liq*(T_out-T_ref) + V*Hvap_mix = cp_liq*(T_in-T_ref) [изоэнтальп.]
    %   => V = cp_mix*(T_in - T_out) / Hvap_mix(T_out)
    % ------------------------------------------------------------------

    if T_in <= T_bub + 0.5
        % Жидкость остаётся жидкостью
        T_out  = T_in;   % изоэнтальп. для несжимаемой жидкости ≈ изотерма
        V_frac = 0;
        phase  = 'liquid';
        y_out  = zeros(1,3);
        x_out  = x;
    else
        % Двухфазная область: решаем итерационно
        T_out = T_bub;   % начальная оценка

        for iter = 1:100
            % Равновесный состав пара
            flash_opts.tol = 1e-8;
            [V_try, y_try, x_try, ~, ~] = vle_flash(x, P_out, T_out, ...
                ant, delta_g, alpha, R, flash_opts);

            if V_try <= 0
                T_out = T_bub; V_frac = 0;
                y_out = zeros(1,3); x_out = x;
                break
            end

            % Теплота испарения при T_out
            Hv = Hvap_mix(T_out, y_try);   % кДж/кмоль

            % Уравнение баланса: cp_mix*(T_in - T_out) = V * Hv
            cp_val = mean([cp_mix(T_in), cp_mix(T_out)]);
            V_new  = cp_val * (T_in - T_out) / max(Hv, 1);

            V_new = min(max(V_new, 0), 1);

            % Уточнение T_out из условия V_new = V_try
            if abs(V_new - V_try) < 1e-5
                V_frac = V_try;
                y_out  = y_try;
                x_out  = x_try;
                break
            end
            % Обновление T_out через явную формулу
            T_out_new = T_in - V_new * Hv / max(cp_val, 1);
            T_out_new = max(min(T_out_new, T_in), T_bub - 5);
            T_out = 0.5*T_out_new + 0.5*T_out;
        end

        if T_out < T_bub
            T_out = T_bub;
        end

        if V_frac > 0.99
            phase = 'vapor';
        else
            phase = 'two-phase';
        end
    end

    % ------------------------------------------------------------------
    % Коэффициент Джоуля-Томсона (приближение)
    % mu_JT = (1/cp) * (T*dV/dT - V) ≈ -ΔT/ΔP [К/Па]
    % ------------------------------------------------------------------
    dP = P_in - P_out;
    dT = T_in - T_out;
    if abs(dP) > 1e3
        mu_JT = dT / dP;
    else
        mu_JT = 0;
    end

    % ------------------------------------------------------------------
    % Формирование выходных потоков
    % ------------------------------------------------------------------
    n_vap = n * V_frac;
    n_liq = n * (1 - V_frac);

    if ~exist('y_out','var'), y_out = zeros(1,3); end
    if ~exist('x_out','var'), x_out = x; end

    result.T_out  = T_out;
    result.P_out  = P_out;
    result.phase  = phase;
    result.V_frac = V_frac;
    result.T_bub  = T_bub;
    result.dT     = T_in - T_out;
    result.Joule_Thomson = mu_JT;

    result.vapor.n  = n_vap;
    result.vapor.x  = y_out;
    result.vapor.T  = T_out;
    result.vapor.P  = P_out;

    result.liquid.n = n_liq;
    result.liquid.x = x_out;
    result.liquid.T = T_out;
    result.liquid.P = P_out;

    name = getf(feed,'name','Клапан');
    fprintf('  %-20s  P_in=%6.1f→%6.1f кПа  T_in=%6.2f°C  T_out=%6.2f°C  ΔT=%.3f К\n', ...
        name, P_in/1e3, P_out/1e3, T_in-273.15, T_out-273.15, T_in-T_out);
    fprintf('    Фаза: %-12s  V_frac=%.4f  n_пар=%.2f  n_жидк=%.2f кмоль/ч\n', ...
        phase, V_frac, n_vap, n_liq);
    fprintf('    μ_ДТ = %.3e К/Па\n', mu_JT);
end

function T_b = bubble_T_simple(z, P, ant, delta_g, alpha, R)
    T_b = 320;
    for k = 1:300
        Ps_A=exp(ant(1,1)+ant(1,2)/(T_b+ant(1,3)));
        Ps_B=exp(ant(2,1)+ant(2,2)/(T_b+ant(2,3)));
        gm = nrtl_gamma(z,delta_g,alpha,R,T_b);
        f  = gm(1)*Ps_A/P*z(1)+gm(2)*Ps_B/P*z(2)-1;
        dT = 0.1;
        Ps_A2=exp(ant(1,1)+ant(1,2)/(T_b+dT+ant(1,3)));
        Ps_B2=exp(ant(2,1)+ant(2,2)/(T_b+dT+ant(2,3)));
        gm2=nrtl_gamma(z,delta_g,alpha,R,T_b+dT);
        f2=gm2(1)*Ps_A2/P*z(1)+gm2(2)*Ps_B2/P*z(2)-1;
        df=(f2-f)/dT; if abs(df)<1e-15, break; end
        T_n=T_b-f/df; if abs(T_n-T_b)<1e-6, T_b=T_n; break; end
        T_b=T_n;
    end
end

function val = getf(s,f,d)
    if isfield(s,f), val=s.(f); else, val=d; end
end
