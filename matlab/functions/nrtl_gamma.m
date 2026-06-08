function gamma = nrtl_gamma(x, delta_g, alpha, R, T)
%NRTL_GAMMA  Коэффициенты активности по модели NRTL (n-компонентная смесь)
%
%  Входные параметры:
%    x       — вектор мольных долей [1×n], сумма = 1
%    delta_g — матрица энергий взаимодействия [кДж/моль], n×n, delta_g(i,i)=0
%    alpha   — параметр нерандомности (скаляр; одинаков для всех пар)
%    R       — газовая постоянная [Дж/моль/К]
%    T       — температура [К]
%
%  Выходные параметры:
%    gamma   — коэффициенты активности [1×n], gamma(i) >= 1 (обычно)
%
%  Уравнения NRTL (Chen, 1982; уравнения 14-17 диссертации):
%    tau_ij  = delta_g_ij * 1000 / (R*T)             — безразмерная энергия
%    G_ij    = exp(-alpha * tau_ij)                   — матрица Больцмана
%    S_j     = sum_k( x_k * G_kj )                   — нормировочная сумма
%    ln_gamma_i = sum_j(x_j*tau_ji*G_ji)/S_i
%               + sum_j[ x_j*G_ij/S_j * (tau_ij - sum_k(x_k*tau_kj*G_kj)/S_j) ]
%
%  Пример:
%    delta_g = [0 5.57 19.27; 7.50 0 10.65; 16.21 0.56 0];
%    gamma = nrtl_gamma([0.54 0.46 0.0], delta_g, 0.4, 8.314, 318.15)

    % --- Проверка входных данных ---
    x = max(x(:)', 1e-14);      % строка, защита от нулей
    x = x / sum(x);             % нормировка
    n = length(x);

    % --- Параметры модели ---
    tau = (delta_g * 1000) ./ (R * T);   % матрица tau_ij [безразм.]
    G   = exp(-alpha .* tau);             % матрица G_ij

    % --- Нормировочные суммы S_j = sum_k(x_k * G_kj) ---
    S  = x * G;                   % вектор [1×n]

    % --- Вспомогательные матрицы ---
    xG  = x .* G;                 % x_k * G_kj
    xtG = x .* (tau .* G);        % x_k * tau_kj * G_kj

    % --- Расчёт ln(gamma_i) ---
    ln_gamma = zeros(1, n);
    for i = 1:n
        % Первое слагаемое: sum_k(x_k*tau_ki*G_ki) / S_i
        term1 = sum(xtG(:,i)') / S(i);

        % Второе слагаемое: sum_j[ x_j*G_ij/S_j * (tau_ij - numj/S_j) ]
        term2 = 0;
        for j = 1:n
            numj  = sum(xtG(:,j)');   % sum_k(x_k*tau_kj*G_kj)
            term2 = term2 + (xG(j) * G(i,j) / S(j)) * (tau(i,j) - numj/S(j));
        end

        ln_gamma(i) = term1 + term2;
    end

    gamma = exp(ln_gamma);
end
