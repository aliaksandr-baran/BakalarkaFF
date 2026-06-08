

%% Initialization
clc
clear all
close all

%% LLE - equlibrium and NRTL model
           % A      B     C
xR = [0.90, 0.05, 0.05]; % Rafinat - nastrel
xE = [0.05, 0.90, 0.05]; % Extrakt - nastrel

tol = 1e-8;
error = 1;
max_iter = 1000;
iter = 0;

% Parameters
R = 8.314; % [J/mol/K]
T = 318.15; % [K]
t = T-273.15; % [C]

alpha = 0.4;
% [11 12 13; 21 22 23; 31 32 33]
delta_g = [0.0000 5.5737 19.2694;
              7.4991 0.0000 0.5632; 
              16.2120 0.5632 0.000]; % [kJ/mol]

% tau = delta_g./ (R/1000*T);
% G = exp(-alpha.*tau);

% Iteračný cyklus (Successive Substitution)
while error > tol && iter < max_iter
    iter = iter + 1;

    % Výpočet aktivitných koeficientov pre obe fázy
    gamma_R = nrtl_gamma(xR, delta_g, alpha, R, T);
    gamma_E  = nrtl_gamma(xE, delta_g, alpha, R, T);
    
    % Výpočet distribučných koeficientov K_i
    K= gamma_R ./ gamma_E;
    % K_R2E = gamma_R ./ gamma_E;
    % K_E2R = gamma_E ./ gamma_R;

    xE_old = xE;
    xE = xE ./ sum(K .* xE);
    xR = K .* xE;
    xR = xR / sum(xR);

    error = sum(abs(xE - xE_old));
end

fprintf('Rovnováha dosiahnutá po %d iteráciách:\n', iter);
fprintf('Fáza alfa: [%.4f, %.4f, %.4f]\n', xR);
fprintf('Fáza beta:  [%.4f, %.4f, %.4f]\n', xE);


% Pomocná funkcia
function gamma = nrtl_gamma(x, delta_g, alpha, R, T)
    n = length(x);
    tau = delta_g.*1000 ./ (R * T);
    G = exp(-alpha .* tau);
    S = x .* G;
    term1 = (x .* (tau .* G)) ./ S;
    term2 = ((x .* G) ./ S) .* (tau - (ones(n,1) .* ((x .* (tau .* G)) ./ S)))';
    gamma = exp(term1 + term2);
end
