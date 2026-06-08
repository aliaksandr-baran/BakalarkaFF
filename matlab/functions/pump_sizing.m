function result = pump_sizing(stream, P_in, P_out, opts)
%PUMP_SIZING  Расчёт насоса: гидравлика, мощность, подбор типа
%
%  Входные параметры:
%    stream — структура потока:
%               .n    — молярный поток [кмоль/ч]
%               .x    — мольные доли [xA, xB, xC]
%               .T    — температура [К]
%               .name — имя (строка)
%    P_in   — давление на всасывании [Па]
%    P_out  — давление на нагнетании [Па]
%    opts   — (необязательно):
%               .eta_pump — КПД насоса (0.70)
%               .eta_motor — КПД двигателя (0.95)
%               .dz       — геодезический подъём [м] (1.0)
%               .g        — ускорение свободного падения (9.81 м/с²)
%               .rho_A, .rho_B, .rho_C — плотности [кг/м³]
%               .MA, .MB, .MC — молярные массы [кг/кмоль]
%
%  Выходные параметры:
%    result — структура:
%      .Q_vol      — объёмный расход [м³/ч]
%      .rho_mix    — плотность смеси [кг/м³]
%      .dP         — перепад давления [Па]
%      .H_head     — напор [м вод.ст.]
%      .P_hydraulic — гидравлическая мощность [Вт]
%      .P_shaft    — мощность на валу [Вт] (с учётом КПД насоса)
%      .P_motor    — мощность двигателя [кВт] (с учётом КПД двигателя)
%      .NPSH_req   — требуемый кавитационный запас [м]
%      .pump_type  — рекомендуемый тип насоса (строка)

    if nargin < 4, opts = struct(); end
    eta_p  = getf(opts,'eta_pump',  0.70);
    eta_m  = getf(opts,'eta_motor', 0.95);
    dz     = getf(opts,'dz',        1.0);
    g      = getf(opts,'g',         9.81);

    % Молярные массы и плотности
    MA = getf(opts,'MA', 88.15);  MB = getf(opts,'MB', 32.04);  MC = getf(opts,'MC', 236.29);
    rho_A = getf(opts,'rho_A', 740);    % кг/м³ (MTBE при 45°C)
    rho_B = getf(opts,'rho_B', 775);    % кг/м³ (MeOH при 45°C)
    rho_C = getf(opts,'rho_C', 1275);   % кг/м³ ([BMIM][HSO4] при 45°C)

    x   = max(stream.x(:)', 0); x = x/sum(x);
    n   = stream.n;   % кмоль/ч

    % --- Молярная масса смеси [кг/кмоль] ---
    M_mix = x(1)*MA + x(2)*MB + x(3)*MC;

    % --- Плотность смеси по правилу смешения объёмов (приближение) ---
    % 1/rho_mix = sum(w_i/rho_i), где w_i — массовые доли
    w_A = x(1)*MA / M_mix;
    w_B = x(2)*MB / M_mix;
    w_C = x(3)*MC / M_mix;
    rho_mix = 1 / (w_A/rho_A + w_B/rho_B + w_C/rho_C);

    % --- Массовый и объёмный расходы ---
    G_mass = n * M_mix;            % кг/ч
    Q_vol  = G_mass / rho_mix;     % м³/ч
    Q_vol_s = Q_vol / 3600;        % м³/с

    % --- Гидравлические параметры ---
    dP = P_out - P_in;             % Па
    H_head = dP / (rho_mix * g) + dz;   % м (напор полный)

    % --- Мощности ---
    P_hydr  = Q_vol_s * dP;               % Вт (гидравлическая)
    P_shaft = P_hydr / eta_p;             % Вт (на валу насоса)
    P_motor = P_shaft / eta_m / 1000;     % кВт (электрическая)

    % --- Требуемый кавитационный запас (NPSH_req) ---
    % Эмпирическая корреляция: NPSH_req = 0.3 + 0.03*Q^0.5 [м]
    NPSH_req = 0.3 + 0.03 * sqrt(Q_vol);

    % --- Рекомендуемый тип насоса ---
    pump_type = select_pump_type(Q_vol, dP);

    % --- Результат ---
    result.Q_vol       = Q_vol;
    result.rho_mix     = rho_mix;
    result.M_mix       = M_mix;
    result.G_mass      = G_mass;
    result.dP          = dP;
    result.H_head      = H_head;
    result.P_hydraulic = P_hydr;
    result.P_shaft     = P_shaft;
    result.P_motor     = P_motor;
    result.NPSH_req    = NPSH_req;
    result.pump_type   = pump_type;
    result.eta_pump    = eta_p;
    result.eta_motor   = eta_m;

    name = getf(stream,'name','Насос');
    fprintf('  %-12s  Q=%6.2f м³/ч  dP=%7.1f кПа  H=%6.1f м  P_вал=%7.1f Вт  P_дв=%5.2f кВт  [%s]\n', ...
        name, Q_vol, dP/1e3, H_head, P_shaft, P_motor, pump_type);
end

% -----------------------------------------------------------------------
function ptype = select_pump_type(Q_vol, dP)
%SELECT_PUMP_TYPE  Рекомендация типа насоса по расходу и напору
    dP_bar = dP / 1e5;
    if Q_vol < 5 && dP_bar > 10
        ptype = 'плунжерный';
    elseif Q_vol < 20
        ptype = 'центробежный (малый)';
    elseif Q_vol < 200
        ptype = 'центробежный';
    else
        ptype = 'центробежный (крупный)';
    end
end

function val = getf(s,f,d)
    if isfield(s,f), val=s.(f); else, val=d; end
end
