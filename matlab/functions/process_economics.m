function result = process_economics(equipment_list, utility_rates, annual_hours, opts)
%PROCESS_ECONOMICS  CAPEX / OPEX расчёт с индексированием CEPCI
%
%  Реализует методологию оценки стоимости оборудования по Guthrie/Peters-Timmerhaus:
%    C_eq(2025) = C_eq(base) * (CEPCI_2025 / CEPCI_base)
%    C_installed = C_eq * f_install (коэффициент монтажа)
%    TCI = C_installed_total * f_TCI (полные капиталовложения)
%
%  Входные параметры:
%    equipment_list — массив структур с полями:
%                       .name     — имя (строка)
%                       .type     — 'pump','HEX','column','vessel','other'
%                       .cost_base — базовая стоимость [€] при CEPCI_base
%                       .CEPCI_base — базовый CEPCI (по умолч. 567.5 = 2017)
%                       .size_param — характерный размер (площадь, расход и т.д.)
%    utility_rates  — структура:
%                       .steam_price  — стоимость пара [€/т]
%                       .cooling_price— стоимость охл. воды [€/т]
%                       .elec_price   — стоимость эл. энергии [€/кВт·ч]
%                       .steam_kg_h   — расход пара [т/ч]
%                       .cooling_kg_h — расход охл. воды [т/ч]
%                       .elec_kW      — потребление эл. энергии [кВт]
%    annual_hours   — число рабочих часов в год
%    opts           — (необязательно):
%                       .CEPCI_now   — текущий CEPCI (811.5 = 2025)
%                       .f_install   — коэффициент монтажа (3.1)
%                       .f_TCI       — коэффициент ТСИ (1.18)
%                       .f_labor     — доля затрат труда от TCI/год (0.02)
%                       .f_maint     — доля затрат на ТО от TCI/год (0.01)
%                       .IL_cost     — стоимость ионной жидкости [€]
%
%  Выходные параметры:
%    result — структура:
%      .equip_costs  — скорректированные стоимости оборудования [€]
%      .CAPEX_equip  — суммарные затраты на оборудование [€]
%      .CAPEX_install— затраты с учётом монтажа [€]
%      .TCI          — полные капиталовложения (TCI) [€]
%      .OPEX         — операционные расходы [€/год]
%      .OPEX_breakdown — разбивка OPEX по статьям
%      .payback_years— простой срок окупаемости (если задана прибыль) [лет]

    if nargin < 4, opts = struct(); end
    CEPCI_now  = getf(opts,'CEPCI_now',  811.5);   % CEPCI 2025 (оценка)
    f_install  = getf(opts,'f_install',  3.1);     % Lang-фактор для химпроизводства
    f_TCI      = getf(opts,'f_TCI',      1.18);    % накладные, разрешения
    f_labor    = getf(opts,'f_labor',    0.02);    % 2% от TCI в год (труд)
    f_maint    = getf(opts,'f_maint',    0.01);    % 1% от TCI в год (ТО)
    IL_cost    = getf(opts,'IL_cost',    0);

    n_eq = length(equipment_list);
    equip_costs = zeros(1, n_eq);

    fprintf('\n--- CAPEX: Стоимость оборудования ---\n');
    fprintf('%-22s  %12s  %12s  %12s\n', 'Оборудование', 'C_base [€]', 'CEPCI-индекс', 'C_2025 [€]');
    fprintf('%s\n', repmat('-',1,64));

    for i = 1:n_eq
        eq = equipment_list(i);
        cb  = eq.cost_base;
        ci_base = getf(eq,'CEPCI_base', 567.5);
        c_now = cb * (CEPCI_now / ci_base);
        equip_costs(i) = c_now;
        fprintf('%-22s  %12.0f  %12.1f  %12.0f\n', eq.name, cb, ci_base, c_now);
    end

    CAPEX_equip   = sum(equip_costs) + IL_cost;
    CAPEX_install = CAPEX_equip * f_install;
    TCI           = CAPEX_install * f_TCI;

    fprintf('%s\n', repmat('-',1,64));
    fprintf('%-22s  %36.0f\n', 'Итого оборудование:', CAPEX_equip);
    fprintf('%-22s  %36.0f\n', sprintf('С монтажом (f=%.1f):',f_install), CAPEX_install);
    fprintf('%-22s  %36.0f\n', 'TCI:', TCI);

    % --- OPEX ---
    cost_steam   = utility_rates.steam_price   * utility_rates.steam_kg_h  * annual_hours / 1000;
    cost_cooling = utility_rates.cooling_price * utility_rates.cooling_kg_h* annual_hours / 1000;
    cost_elec    = utility_rates.elec_price    * utility_rates.elec_kW     * annual_hours;
    cost_labor   = TCI * f_labor;
    cost_maint   = TCI * f_maint;
    OPEX_total   = cost_steam + cost_cooling + cost_elec + cost_labor + cost_maint;

    fprintf('\n--- OPEX: Операционные расходы [€/год] ---\n');
    fprintf('  Греющий пар:         %12.0f\n', cost_steam);
    fprintf('  Охлаждающая вода:    %12.0f\n', cost_cooling);
    fprintf('  Электроэнергия:      %12.0f\n', cost_elec);
    fprintf('  Труд персонала:      %12.0f\n', cost_labor);
    fprintf('  Техническое обслуж.: %12.0f\n', cost_maint);
    fprintf('  ИТОГО OPEX:          %12.0f  €/год\n', OPEX_total);

    % --- Результат ---
    result.equip_costs    = equip_costs;
    result.CAPEX_equip    = CAPEX_equip;
    result.CAPEX_install  = CAPEX_install;
    result.TCI            = TCI;
    result.OPEX           = OPEX_total;
    result.OPEX_breakdown = struct(...
        'steam',   cost_steam, ...
        'cooling', cost_cooling, ...
        'elec',    cost_elec, ...
        'labor',   cost_labor, ...
        'maint',   cost_maint);

    % --- График ---
    figure('Name','Экономика процесса','Position',[50 50 900 400]);

    subplot(1,2,1);
    names_plot = {equipment_list.name};
    vals_plot  = equip_costs / 1e3;
    [sv, si] = sort(vals_plot,'descend');
    bar(sv,'FaceColor',[0.3 0.6 0.9]);
    set(gca,'XTick',1:n_eq,'XTickLabel',names_plot(si));
    xtickangle(40); ylabel('тыс. €');
    title(sprintf('CAPEX по статьям\nИТОГО: %.0f тыс. €', CAPEX_equip/1e3));
    grid on;

    subplot(1,2,2);
    opex_labels = {'Пар','Охл. вода','Эл.энергия','Труд','ТО'};
    opex_vals   = [cost_steam, cost_cooling, cost_elec, cost_labor, cost_maint];
    pie(opex_vals);
    legend(opex_labels,'Location','south','FontSize',8);
    title(sprintf('OPEX по статьям\nИТОГО: %.0f тыс. €/год', OPEX_total/1e3));
end

function val = getf(s,f,d)
    if isfield(s,f), val=s.(f); else, val=d; end
end
