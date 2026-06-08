function result = mass_transfer_analysis(stage_xR, stage_xE, stage_nR, stage_nE, ...
                                          delta_g, alpha, R, T, opts)
%MASS_TRANSFER_ANALYSIS  Детальный анализ массопереноса в экстракторе
%
%  Рассчитывает для каждой ступени:
%    - коэффициенты распределения K_i
%    - движущие силы (∆x_i, ∆y_i) — разность фактического и равновесного составов
%    - эффективность Мёрфри EMV_i по каждому компоненту
%    - локальные числа NTU и HTU
%    - вклад каждой ступени в общий массоперенос
%
%  Входные параметры:
%    stage_xR  — профиль рафинатов  [(N+1)×3]
%    stage_xE  — профиль экстрактов [(N+1)×3]
%    stage_nR  — молярные потоки рафинатов [кмоль/ч], вектор (N+1)
%    stage_nE  — молярные потоки экстрактов [кмоль/ч], вектор (N+1)
%    delta_g   — матрица NRTL [кДж/моль]
%    alpha, R, T — параметры NRTL
%    opts      — (необязательно):
%                  .HETS  — высота теор. ступени [м] (0.15)
%                  .comp_names — {'MTBE','MeOH','IL'}
%
%  Выходные параметры:
%    result — структура:
%      .K             — коэффициенты распределения [N×3]
%      .driving_force_R — движ. сила по рафинату [N×3]
%      .driving_force_E — движ. сила по экстракту [N×3]
%      .E_Murphree_R  — эффективность Мёрфри по рафинату [N×3]
%      .E_Murphree_E  — эффективность Мёрфри по экстракту [N×3]
%      .lambda        — фактор экстракции (stripping factor) [N×3]
%      .NTU_per_stage — число ед. переноса на ступень [N×1]
%      .HTU_per_stage — высота ед. переноса [м] [N×1]
%      .NTU_total     — суммарное NTU
%      .selectivity   — селективность β_BA по ступеням [N×1]

    if nargin < 9, opts = struct(); end
    HETS = getf(opts,'HETS', 0.15);   % м
    comp_names = getf(opts,'comp_names', {'MTBE','MeOH','IL'});

    N = size(stage_xR,1) - 1;   % число ступеней

    K_mat    = zeros(N, 3);
    DF_R     = zeros(N, 3);   % движущая сила по рафинату (x - x*)
    DF_E     = zeros(N, 3);   % движущая сила по экстракту (y* - y)
    E_MV_R   = zeros(N, 3);
    E_MV_E   = zeros(N, 3);
    lambda   = zeros(N, 3);   % фактор экстракции = K * L/V
    NTU_stage= zeros(N, 1);
    sel      = zeros(N, 1);

    fprintf('\n--- Анализ массопереноса по ступеням экстрактора ---\n');
    fprintf('%-5s  %-6s  %-8s  %-8s  %-8s  %-8s  %-8s  %-8s  %-6s  %-5s\n', ...
        'Ступ', 'K_A', 'K_B', 'K_C', 'ΔxA', 'ΔxB', 'EMV_A', 'EMV_B', 'λ_B', 'β_BA');
    fprintf('%s\n', repmat('-',1,82));

    for s = 1:N
        xR_in  = stage_xR(s,:);    % рафинат, входящий в ступень
        xR_out = stage_xR(s+1,:);  % рафинат, выходящий из ступени
        xE_in  = stage_xE(s+1,:);  % экстракт, входящий в ступень (снизу)
        xE_out = stage_xE(s,:);    % экстракт, выходящий из ступени (сверху)

        % Коэффициенты активности и K_i = gamma_E / gamma_R
        gR = nrtl_gamma(xR_out, delta_g, alpha, R, T);
        gE = nrtl_gamma(xE_out, delta_g, alpha, R, T);
        K  = gE ./ gR;
        K_mat(s,:) = K;

        % --- Равновесный состав рафината (при данном экстракте) ---
        xR_eq = K .* xE_out;
        if sum(xR_eq) > 0, xR_eq = xR_eq/sum(xR_eq); end

        % --- Равновесный состав экстракта (при данном рафинате) ---
        xE_eq = xR_out ./ max(K, 1e-15);
        if sum(xE_eq) > 0, xE_eq = xE_eq/sum(xE_eq); end

        % --- Движущие силы ---
        DF_R(s,:) = xR_out - xR_eq;    % x - x* (должна быть > 0 если MTBE извлекается)
        DF_E(s,:) = xE_eq  - xE_out;   % y* - y

        % --- Эффективность Мёрфри по рафинату ---
        % EMV = (x_in - x_out) / (x_in - x*_out)
        denom_R = xR_in - xR_eq;
        for ci = 1:3
            if abs(denom_R(ci)) > 1e-8
                E_MV_R(s,ci) = (xR_in(ci) - xR_out(ci)) / denom_R(ci);
            else
                E_MV_R(s,ci) = 1.0;
            end
        end

        % --- Эффективность Мёрфри по экстракту ---
        denom_E = xE_eq - xE_in;
        for ci = 1:3
            if abs(denom_E(ci)) > 1e-8
                E_MV_E(s,ci) = (xE_out(ci) - xE_in(ci)) / denom_E(ci);
            else
                E_MV_E(s,ci) = 1.0;
            end
        end

        % --- Фактор экстракции λ_i = K_i * (L_E / L_R) ---
        L_R = max(stage_nR(s), 1e-9);
        L_E = max(stage_nE(s+1), 1e-9);
        lambda(s,:) = K * (L_E / L_R);

        % --- Локальное NTU (по компоненту B = MeOH, ключевому) ---
        avg_df = max(abs(DF_R(s,2)), 1e-9);
        NTU_stage(s) = abs(xR_in(2) - xR_out(2)) / avg_df;

        % --- Селективность β_BA ---
        K_B = K(2); K_A = K(1);
        if abs(K_A) > 1e-9
            sel(s) = K_B / K_A;
        else
            sel(s) = NaN;
        end

        fprintf('%-5d  %-6.3f  %-8.3f  %-8.4f  %-8.5f  %-8.5f  %-8.3f  %-8.3f  %-6.3f  %-5.2f\n', ...
            s, K(1), K(2), K(3), DF_R(s,1), DF_R(s,2), E_MV_R(s,1), E_MV_R(s,2), lambda(s,2), sel(s));
    end

    NTU_total = sum(NTU_stage);
    HTU_stage = HETS * ones(N,1) ./ max(NTU_stage, 1e-6);   % HTU = HETS/NTU_local
    HTU_stage(NTU_stage < 1e-6) = NaN;

    fprintf('%s\n', repmat('-',1,82));
    fprintf('  Суммарное NTU (по MeOH): %.3f\n', NTU_total);
    fprintf('  Средняя эффективность Мёрфри EMV_MeOH: %.1f%%\n', ...
        mean(E_MV_R(isfinite(E_MV_R(:,2)),2))*100);
    fprintf('  Средняя селективность β_BA: %.2f\n', nanmean(sel));

    result.K              = K_mat;
    result.driving_force_R= DF_R;
    result.driving_force_E= DF_E;
    result.E_Murphree_R   = E_MV_R;
    result.E_Murphree_E   = E_MV_E;
    result.lambda         = lambda;
    result.NTU_per_stage  = NTU_stage;
    result.HTU_per_stage  = HTU_stage;
    result.NTU_total      = NTU_total;
    result.selectivity    = sel;

    % --- Графики ---
    plot_mass_transfer(N, K_mat, E_MV_R, DF_R, lambda, sel);
end

% -----------------------------------------------------------------------
function plot_mass_transfer(N, K, EMV, DF, lambda, sel)
    stages = 1:N;
    figure('Name','Анализ массопереноса','Position',[100 100 1100 750]);

    subplot(2,3,1);
    semilogy(stages, K(:,1),'bo-', stages, K(:,2),'rs-', stages, K(:,3),'g^-','LineWidth',2);
    xlabel('Ступень'); ylabel('K_i (лог. шкала)');
    title('Коэффициенты распределения K_i'); grid on;
    legend({'K_{MTBE}','K_{MeOH}','K_{IL}'},'Location','best');
    xticks(stages);

    subplot(2,3,2);
    bar(stages, EMV(:,1:2)*100, 'grouped');
    xlabel('Ступень'); ylabel('EMV [%]');
    title('Эффективность Мёрфри EMV'); grid on;
    legend({'MTBE','MeOH'},'Location','best'); ylim([0 120]);
    xticks(stages);

    subplot(2,3,3);
    plot(stages, DF(:,1),'bo-', stages, DF(:,2),'rs-','LineWidth',2);
    xlabel('Ступень'); ylabel('Движ. сила Δx_i');
    title('Движущая сила массопереноса (рафинат)'); grid on;
    legend({'MTBE','MeOH'},'Location','best');
    xticks(stages);

    subplot(2,3,4);
    bar(stages, lambda(:,1:2), 'grouped');
    xlabel('Ступень'); ylabel('λ_i = K_i·(L_E/L_R)');
    title('Фактор экстракции λ_i'); grid on;
    legend({'λ_{MTBE}','λ_{MeOH}'},'Location','best');
    xticks(stages);

    subplot(2,3,5);
    plot(stages, sel,'mo-','LineWidth',2,'MarkerSize',8);
    xlabel('Ступень'); ylabel('β_{BA} = K_B/K_A');
    title('Селективность β_{MeOH/MTBE} по ступеням'); grid on;
    xticks(stages);

    subplot(2,3,6);
    % Состав экстракта и рафината по ступеням на одном графике
    text(0.5,0.5,{'Профиль использует','данные из extractor\_column','(см. График 2)'},...
        'HorizontalAlignment','center','FontSize',10,'Units','normalized');
    axis off;
    title('Справка');

    sgtitle('Детальный анализ массопереноса в экстракционной колонне', ...
        'FontSize',13,'FontWeight','bold');
end

function val = getf(s,f,d)
    if isfield(s,f), val=s.(f); else, val=d; end
end
