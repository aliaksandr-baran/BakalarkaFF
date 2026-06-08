function [xE, xR, beta, iter, flag] = lle_solver(z, delta_g, alpha, R, T, opts)
%LLE_SOLVER  Расчёт жидкость-жидкость равновесия (ЖЖЭ)
%            методом последовательной подстановки с демпфированием
%
%  Входные параметры:
%    z       — общий состав системы [1×3] (мольные доли A, B, C)
%    delta_g — матрица NRTL [кДж/моль], 3×3
%    alpha   — параметр нерандомности NRTL
%    R       — газовая постоянная [Дж/моль/К]
%    T       — температура [К]
%    opts    — (необязательно) структура с полями:
%                .tol      — допуск сходимости (по умолчанию 1e-10)
%                .max_iter — макс. итераций (по умолчанию 2000)
%                .damp     — коэффициент демпфирования [0,1] (по умолчанию 0.2)
%                .xR0, .xE0 — начальные составы фаз (по умолчанию [0.97,0.02,0.01])
%
%  Выходные параметры:
%    xE   — мольные доли в экстрактной  фазе (богатой IL) [1×3]
%    xR   — мольные доли в рафинатной   фазе (богатой MTBE) [1×3]
%    beta — доля экстрактной фазы (мольная) [безразм.]
%    iter — число выполненных итераций
%    flag — 0: сошлось и нетривиально; 1: не сошлось; 2: тривиальное решение
%
%  Метод:
%    Уравнение равновесия: gamma_i^R * x_i^R = gamma_i^E * x_i^E
%    => K_i = x_i^R / x_i^E = gamma_i^E / gamma_i^R
%    Последовательная замена с демпфированием для устойчивости:
%      xE_new = xE / sum(K .* xE)
%      xE = damp*xE_new + (1-damp)*xE_old
%
%  Пример:
%    [xE, xR, beta] = lle_solver([0.3 0.3 0.4], delta_g, 0.4, 8.314, 318.15)

    % --- Параметры по умолчанию ---
    if nargin < 6 || isempty(opts), opts = struct(); end
    tol      = getfield_def(opts, 'tol',      1e-10);
    max_iter = getfield_def(opts, 'max_iter', 2000);
    damp     = getfield_def(opts, 'damp',     0.2);
    xR0      = getfield_def(opts, 'xR0',      [0.97, 0.02, 0.01]);
    xE0      = getfield_def(opts, 'xE0',      [0.01, 0.02, 0.97]);

    xR = xR0 / sum(xR0);
    xE = xE0 / sum(xE0);

    diff_val = 1;
    iter     = 0;
    flag     = 1;

    % --- Итерационный цикл ---
    while diff_val > tol && iter < max_iter
        iter = iter + 1;

        gR = nrtl_gamma(xR, delta_g, alpha, R, T);
        gE = nrtl_gamma(xE, delta_g, alpha, R, T);

        % Коэффициенты распределения K_i = gamma_E_i / gamma_R_i
        K = gE ./ gR;

        xR_old = xR;
        xE_old = xE;

        % Обновление экстрактной фазы
        xE_new = xE ./ sum(K .* xE);
        xE_new = max(xE_new, 0);
        xE_new = xE_new / sum(xE_new);

        xE = damp * xE_new + (1 - damp) * xE_old;
        xE = max(xE, 0); xE = xE / sum(xE);

        % Рафинатная фаза через K
        xR = K .* xE;
        xR = max(xR, 0); xR = xR / sum(xR);

        diff_val = max(abs(xR - xR_old)) + max(abs(xE - xE_old));
    end

    % --- Доля экстрактной фазы (уравнение Ричфорда-Райса) ---
    beta = lle_phase_fraction(z, xR, xE);

    % --- Классификация результата ---
    separation = sum(abs(xR - xE));
    if diff_val <= tol
        if separation > 5e-3
            flag = 0;   % сошлось, нетривиально
        else
            flag = 2;   % тривиальное решение (гомогенная фаза)
        end
    end
end

% -----------------------------------------------------------------------
function beta = lle_phase_fraction(z, xR, xE)
% Доля экстрактной фазы из баланса: z = beta*xE + (1-beta)*xR
% Решение по МНК (один компонент даёт точный ответ, усредняем по всем)
    d = xE - xR;
    mask = abs(d) > 1e-6;
    if any(mask)
        beta = mean((z(mask) - xR(mask)) ./ d(mask));
        beta = min(max(beta, 0), 1);
    else
        beta = 0.5;
    end
end

% -----------------------------------------------------------------------
function val = getfield_def(s, field, default)
    if isfield(s, field)
        val = s.(field);
    else
        val = default;
    end
end
