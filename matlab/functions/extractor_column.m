function result = extractor_column(feed, solvent, N_stages, delta_g, alpha, R, T, opts)
%EXTRACTOR_COLUMN  Противоточный жидкость-жидкость экстрактор
%                  Материальный баланс + масса- и теплоперенос
%
%  Входные параметры:
%    feed    — структура питания:
%                .n  — молярный поток [кмоль/ч]
%                .x  — мольные доли [xA, xB, xC]
%                .T  — температура [К]
%    solvent — структура растворителя (те же поля)
%    N_stages — число теоретических ступеней
%    delta_g  — матрица NRTL [кДж/моль]
%    alpha    — параметр нерандомности NRTL
%    R        — газовая постоянная [Дж/моль/К]
%    T        — рабочая температура экстрактора [К]
%    opts     — (необязательно) структура:
%                .tol      — допуск (1e-8)
%                .max_iter — макс. итераций (3000)
%                .damp     — демпфирование состава (0.15)
%                .E_stage  — эффективность Мёрфри (0.85)
%                .HETS     — высота, эквивалентная теор. ступени [м] (0.15)
%                .d_col    — диаметр колонны [м] (расчёт из скорости)
%                .rho_R    — плотность рафинатной фазы [кг/м³] (750)
%                .rho_E    — плотность экстрактной  фазы [кг/м³] (1230)
%                .sigma    — межфазное натяжение [мН/м] (15)
%
%  Выходные параметры:
%    result — структура:
%      .raffinate    — поток рафината: .n, .x
%      .extract      — поток экстракта: .n, .x
%      .stage_xR     — профиль рафинатов  [(N+1)×3]
%      .stage_xE     — профиль экстрактов [(N+1)×3]
%      .stage_nR     — молярные потоки рафинатов [кмоль/ч]
%      .stage_nE     — молярные потоки экстрактов [кмоль/ч]
%      .K_profile    — коэффициенты распределения по ступеням [N×3]
%      .N_stages_req — число реальных ступеней (= N/E_stage)
%      .col_height   — высота колонны [м]
%      .col_diameter — диаметр колонны [м]
%      .flooding_vel — скорость захлёбывания [м/с]
%      .oper_vel     — рабочая скорость дисперс. фазы [м/с]
%      .iter         — число итераций
%      .error        — остаточная невязка

    % --- Параметры по умолчанию ---
    if nargin < 8, opts = struct(); end
    tol      = getf(opts, 'tol',      1e-8);
    max_iter = getf(opts, 'max_iter', 3000);
    damp     = getf(opts, 'damp',     0.15);
    E_stage  = getf(opts, 'E_stage',  0.85);   % эффективность Мёрфри
    HETS     = getf(opts, 'HETS',     0.15);   % [м] — типично для пульс. колонн с ИЖ
    rho_R    = getf(opts, 'rho_R',    750);    % кг/м³ рафинат (MTBE-богатый)
    rho_E    = getf(opts, 'rho_E',    1230);   % кг/м³ экстракт (ИЖ-богатый)
    sigma    = getf(opts, 'sigma',    15);     % мН/м межфазное натяжение

    N = N_stages;
    n_F = feed.n;    x_F = feed.x;
    n_S = solvent.n; x_S = solvent.x;

    % --- Начальные матрицы ---
    % Нумерация: строка 1 = питание/выход экстракта, строка N+1 = растворитель/рафинат
    xR = zeros(N+1, 3);  xE = zeros(N+1, 3);
    nR = zeros(1, N+1);  nE = zeros(1, N+1);

    xR(1,:) = x_F;    nR(1) = n_F;
    xE(N+1,:) = x_S;  nE(N+1) = n_S;

    % Линейная интерполяция начального профиля
    for i = 2:N+1
        t = (i-1)/N;
        xR(i,:) = (1-t)*x_F + t*[0.99, 0.005, 0.005];
        xR(i,:) = xR(i,:)/sum(xR(i,:));
    end
    for i = 1:N
        t = (N+1-i)/N;
        xE(i,:) = (1-t)*x_S + t*[0.04, 0.55, 0.41];
        xE(i,:) = xE(i,:)/sum(xE(i,:));
    end

    % Начальные потоки из суммарного баланса
    n_total = n_F + n_S;
    x_total = (n_F*x_F + n_S*x_S) / n_total;
    % Оценка n_E(1) из компонентного баланса по компоненту C (ИЖ)
    nE(1) = n_total * x_total(3) / max(xE(1,3), 1e-6);
    nE(1) = min(nE(1), n_total * 0.9);
    nR(N+1) = n_total - nE(1);
    for i = 2:N
        nR(i) = nR(1) + (nR(N+1)-nR(1))*(i-1)/N;
        nE(i) = nE(1) + (nE(N+1)-nE(1))*(N+1-i)/N;
    end

    % ---------------------------------------------------------------
    % ИТЕРАЦИОННЫЙ ЦИКЛ
    % ---------------------------------------------------------------
    err = 1; n_iter = 0;
    K_profile = zeros(N, 3);

    while err > tol && n_iter < max_iter
        n_iter = n_iter + 1;
        xR_old = xR; xE_old = xE;

        for stage = 1:N
            % Индексы: рафинат входит в ступень из позиции stage,
            %          экстракт входит из позиции stage+1
            gR = nrtl_gamma(xR(stage,:),   delta_g, alpha, R, T);
            gE = nrtl_gamma(xE(stage+1,:), delta_g, alpha, R, T);

            % Коэффициент распределения K_i = x_R_i / x_E_i = gamma_E_i / gamma_R_i
            K = gE ./ gR;
            K_profile(stage,:) = K;

            % Состав рафината на выходе ступени (из уравнения равновесия и баланса)
            % xR(stage+1) = K .* xE(stage) (идеальная ступень)
            xR_eq = K .* xE(stage+1,:);
            if sum(xR_eq) > 0, xR_eq = xR_eq / sum(xR_eq); end

            % Поправка на эффективность Мёрфри:
            % xR_real = xR_in + E * (xR_eq - xR_in)
            xR_new = xR(stage,:) + E_stage * (xR_eq - xR(stage,:));
            xR_new = max(xR_new, 0);
            if sum(xR_new) > 0, xR_new = xR_new / sum(xR_new); end

            xR(stage+1,:) = damp*xR_new + (1-damp)*xR(stage+1,:);
            xR(stage+1,:) = max(xR(stage+1,:),0);
            xR(stage+1,:) = xR(stage+1,:)/sum(xR(stage+1,:));

            % Состав экстракта из равновесия с рафинатом
            gR2 = nrtl_gamma(xR(stage+1,:), delta_g, alpha, R, T);
            gE2 = nrtl_gamma(xE(stage,:),   delta_g, alpha, R, T);
            K2  = gR2 ./ gE2;

            xE_eq = xR(stage+1,:) ./ K2;
            if sum(xE_eq) > 0, xE_eq = xE_eq / sum(xE_eq); end

            xE_new = xE(stage+1,:) + E_stage*(xE_eq - xE(stage+1,:));
            xE_new = max(xE_new, 0);
            if sum(xE_new) > 0, xE_new = xE_new / sum(xE_new); end

            xE(stage,:) = damp*xE_new + (1-damp)*xE(stage,:);
            xE(stage,:) = max(xE(stage,:),0);
            xE(stage,:) = xE(stage,:)/sum(xE(stage,:));
        end

        % Пересчёт потоков из общего баланса
        mat = [1, 1; xE(1,1), xR(N+1,1)];
        rhs = [n_F+n_S; n_F*x_F(1)+n_S*x_S(1)];
        if abs(det(mat)) > 1e-10
            sol = mat \ rhs;
            nE(1)   = max(sol(1), 0);
            nR(N+1) = max(sol(2), 0);
        end
        for i = 1:N
            nE(i+1) = nE(1) + nR(1) - nR(i+1);
            nE(i+1) = max(nE(i+1), 0);
        end

        err = sum(abs(xR(:)-xR_old(:))) + sum(abs(xE(:)-xE_old(:)));
    end

    % ---------------------------------------------------------------
    % РАСЧЁТ МАССОПЕРЕНОСА И РАЗМЕРОВ КОЛОННЫ
    % ---------------------------------------------------------------

    % --- Число реальных ступеней ---
    N_real = ceil(N / E_stage);

    % --- Высота колонны ---
    col_height = N_real * HETS;

    % --- Диаметр колонны из скорости захлёбывания ---
    % Корреляция Miyauchi-Oya для пульсационных колонн:
    %   u_flood [м/с] = C * (sigma/30)^0.2 * ((rho_R-rho_E)/rho_E)^0.5
    % где C ≈ 0.035 м/с для насадки Рашига 6 мм
    C_flood = 0.035;
    u_flood = C_flood * (sigma/30)^0.2 * sqrt((rho_R - rho_E)/rho_E);
    u_oper  = 0.4 * u_flood;   % рабочая скорость = 40% от захлёбывания

    % Объёмный поток рафинатной фазы (питание) — в м³/с
    % Молярный поток n_F [кмоль/ч], средняя мол. масса рафинатной фазы
    M_mix_R = x_F(1)*88.15 + x_F(2)*32.04 + x_F(3)*236.29;
    Q_vol_R = n_F * M_mix_R / rho_R / 3600;   % м³/с
    A_col   = Q_vol_R / max(u_oper, 1e-6);
    d_col   = 2 * sqrt(A_col / pi);
    d_col   = max(d_col, 0.05);   % минимум 5 см

    % --- Объёмный коэффициент массопереноса K_La ---
    % Для пульсационной колонны с ИЖ (оценка):
    %   K_La [1/с] = 0.2 * D_mol^0.5 * (u_oper/u_flood)^0.7 / d_drop^2
    % Упрощённая оценка: K_La ~ 0.5e-3 [1/с] (порядок для IL-систем)
    D_mol  = 1.5e-9;   % коэфф. молекулярной диффузии MeOH в MTBE [м²/с]
    d_drop = 2e-3;     % средний диаметр капли [м]
    ratio  = u_oper / max(u_flood, 1e-6);
    K_La   = 0.2 * sqrt(D_mol) * ratio^0.7 / d_drop^2;

    % --- Число единиц переноса NTU ---
    % NTU = K_La * col_height / u_oper
    NTU = K_La * col_height / max(u_oper, 1e-6);

    % --- HTU = col_height / NTU ---
    HTU = col_height / max(NTU, 1e-6);

    % ---------------------------------------------------------------
    % ФОРМИРОВАНИЕ РЕЗУЛЬТАТА
    % ---------------------------------------------------------------
    result.raffinate.n = nR(N+1);
    result.raffinate.x = xR(N+1,:);
    result.extract.n   = nE(1);
    result.extract.x   = xE(1,:);

    result.stage_xR     = xR;
    result.stage_xE     = xE;
    result.stage_nR     = nR;
    result.stage_nE     = nE;
    result.K_profile    = K_profile;

    result.N_stages_req  = N_real;
    result.E_stage       = E_stage;
    result.col_height    = col_height;
    result.col_diameter  = d_col;
    result.flooding_vel  = u_flood;
    result.oper_vel      = u_oper;
    result.K_La          = K_La;
    result.NTU           = NTU;
    result.HTU           = HTU;
    result.HETS          = HETS;

    result.iter  = n_iter;
    result.error = err;

    % --- Печать итогов ---
    fprintf('\n--- Экстракционная колонна (%d теор. ступеней) ---\n', N);
    fprintf('  Итерации: %d  |  ошибка: %.2e\n', n_iter, err);
    fprintf('  Рафинат:  n=%.2f кмоль/ч  xA=%.4f  xB=%.5f  xC=%.5f\n', ...
        nR(N+1), xR(N+1,1), xR(N+1,2), xR(N+1,3));
    fprintf('  Экстракт: n=%.2f кмоль/ч  xA=%.4f  xB=%.4f  xC=%.4f\n', ...
        nE(1), xE(1,1), xE(1,2), xE(1,3));
    fprintf('  Эффективность Мёрфри: %.0f%%  =>  реальных ступеней: %d\n', E_stage*100, N_real);
    fprintf('  Высота колонны:  %.2f м   Диаметр: %.3f м\n', col_height, d_col);
    fprintf('  Скорость захлёб.: %.4f м/с   Рабочая: %.4f м/с\n', u_flood, u_oper);
    fprintf('  K_La = %.2e 1/с   NTU = %.2f   HTU = %.4f м\n', K_La, NTU, HTU);
end

function val = getf(s, f, d)
    if isfield(s,f), val = s.(f); else, val = d; end
end
