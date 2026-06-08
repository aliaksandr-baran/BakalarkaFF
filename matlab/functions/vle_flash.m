function [V_frac, y, x, T_bub, converged] = vle_flash(z, P, T, ant, delta_g, alpha, R, opts)
%VLE_FLASH  Изотермический ЖПЭ-расчёт (равновесие пар-жидкость)
%           Уравнение Рэтчфорда-Райса + закон Рауля с NRTL-поправкой
%
%  Для двухкомпонентной летучей части системы MTBE(A) + MeOH(B)
%  (компонент C = ионная жидкость считается нелетучим)
%
%  Входные параметры:
%    z       — общий состав [zA, zB, zC] (мольные доли)
%    P       — давление [Па]
%    T       — температура [К]
%    ant     — матрица параметров Антуана [2×3]: строки для A и B
%              ln(Psat[Па]) = ant(i,1) + ant(i,2)/(T + ant(i,3))
%    delta_g — матрица NRTL [кДж/моль], 3×3
%    alpha   — параметр нерандомности
%    R       — газовая постоянная [Дж/моль/К]
%    opts    — (необязательно) структура:
%                .tol      — допуск (по умолчанию 1e-9)
%                .max_iter — макс. итераций (по умолчанию 500)
%
%  Выходные параметры:
%    V_frac    — мольная доля пара (0 = только жидкость, 1 = только пар)
%    y         — состав пара [yA, yB, 0] (ИЖ нелетуча)
%    x         — состав жидкости [xA, xB, xC]
%    T_bub     — температура пузырькования при данном P [К] (рассчитывается)
%    converged — true если итерации сошлись
%
%  Метод:
%    K_i = gamma_i * Psat_i / P    (обобщённый закон Рауля)
%    Уравнение Рэтчфорда-Райса: sum_i [ z_i*(K_i-1)/(1+V*(K_i-1)) ] = 0
%    Решается методом Брента; при V<0 — однофазная жидкость, V>1 — пар.
%
%  Пример:
%    ant = [20.84054 -2624.525 -46.15171;   % MTBE
%           23.5347  -3661.468 -32.77001];  % MeOH
%    [V, y, x] = vle_flash([0.04 0.50 0.46], 2450, 318.15, ant, dg, 0.4, 8.314)

    if nargin < 8 || isempty(opts), opts = struct(); end
    tol      = getfield_def(opts, 'tol',      1e-9);
    max_iter = getfield_def(opts, 'max_iter', 500);

    z = max(z(:)', 1e-15); z = z/sum(z);
    converged = false;

    % --- Давления насыщения при T ---
    Psat_A = exp(ant(1,1) + ant(1,2)/(T + ant(1,3)));
    Psat_B = exp(ant(2,1) + ant(2,2)/(T + ant(2,3)));
    Psat   = [Psat_A, Psat_B, 0];   % ИЖ нелетуча

    % --- Начальные оценки K_i (идеальный газ, жидкость с гамма=1) ---
    x0 = z;
    x0(3) = max(z(3), 1e-10);
    gamma0 = nrtl_gamma(x0, delta_g, alpha, R, T);
    K = gamma0 .* Psat / P;
    K(3) = 0;   % ИЖ остаётся в жидкой фазе

    % --- Проверка на однофазность ---
    sum_Kz = sum(K .* z);   % > 1: есть пар
    sum_z_K= sum(z ./ max(K, 1e-15));  % > 1: есть жидкость (без ИЖ)
    if sum_Kz < 1.0
        % Полностью жидкая фаза
        V_frac = 0; y = zeros(1,3); x = z;
        T_bub = bubble_T(z, P, ant, delta_g, alpha, R, tol, max_iter);
        converged = true; return
    end

    % --- Итерации по уравнению Рэтчфорда-Райса ---
    V_lo = 0.0; V_hi = 1.0;
    V = 0.5;

    for iter = 1:max_iter
        % Состав жидкости
        denom = 1 + V * (K - 1);
        denom(denom < 1e-15) = 1e-15;
        x_it  = z ./ denom;
        x_it(3) = max(x_it(3), 0);
        if sum(x_it) > 0, x_it = x_it / sum(x_it); end

        % Пересчёт K с NRTL
        gamma_it = nrtl_gamma(x_it, delta_g, alpha, R, T);
        K_new  = gamma_it .* Psat / P;
        K_new(3) = 0;

        % Уравнение Рэтчфорда-Райса
        RR = sum(z .* (K_new - 1) ./ (1 + V * (K_new - 1)));

        if abs(RR) < tol
            converged = true;
            K = K_new; x_it_final = x_it;
            break
        end

        % Обновление V методом бисекции
        if RR > 0
            V_lo = V;
        else
            V_hi = V;
        end
        V = 0.5 * (V_lo + V_hi);
        K = K_new;
        x_it_final = x_it;
    end

    % --- Финальный состав ---
    if ~exist('x_it_final','var'), x_it_final = x_it; end
    x = x_it_final;
    x = max(x, 0); x = x / sum(x);

    denom = 1 + V * (K - 1);
    denom(denom < 1e-15) = 1e-15;
    y = K .* x;
    y(3) = 0;
    if sum(y) > 0, y = y / sum(y); end

    V_frac = min(max(V, 0), 1);

    % --- Температура пузырькования ---
    T_bub = bubble_T(z, P, ant, delta_g, alpha, R, tol, max_iter);
end

% -----------------------------------------------------------------------
function T_b = bubble_T(z, P, ant, delta_g, alpha, R, tol, max_iter)
%BUBBLE_T  Температура пузырькования методом Ньютона-Рафсона
%  Условие: sum_i(K_i * z_i) = 1   при V=0, x=z

    T_b = 320;   % начальная оценка [К]
    for k = 1:max_iter
        Ps_A = exp(ant(1,1) + ant(1,2)/(T_b + ant(1,3)));
        Ps_B = exp(ant(2,1) + ant(2,2)/(T_b + ant(2,3)));

        gamma = nrtl_gamma(z, delta_g, alpha, R, T_b);
        K_A = gamma(1) * Ps_A / P;
        K_B = gamma(2) * Ps_B / P;

        f = K_A*z(1) + K_B*z(2) - 1;

        % Числовая производная по T
        dT = 0.1;
        Ps_A2 = exp(ant(1,1) + ant(1,2)/(T_b+dT + ant(1,3)));
        Ps_B2 = exp(ant(2,1) + ant(2,2)/(T_b+dT + ant(2,3)));
        g2 = nrtl_gamma(z, delta_g, alpha, R, T_b+dT);
        f2 = g2(1)*Ps_A2/P*z(1) + g2(2)*Ps_B2/P*z(2) - 1;
        df = (f2 - f) / dT;

        if abs(df) < 1e-15, break; end
        T_new = T_b - f/df;
        if abs(T_new - T_b) < tol, T_b = T_new; break; end
        T_b = T_new;
    end
end

function val = getfield_def(s, field, default)
    if isfield(s, field), val = s.(field); else, val = default; end
end
