function [xE, xR, success] = LLE_solver(z_feed, delta_g, alpha, T)
% LLE_SOLVER Vypočíta rovnovážne zloženie dvoch fáz pre daný nastrel (feed)
% z_feed: Celkové zloženie v bode [z1, z2, z3]
% delta_g: Matica interakčných parametrov [kJ/mol]
% alpha: Parameter neusporiadanosti
% T: Teplota [K]

    R = 8.314;
    tol = 1e-9;
    max_iter = 1500;
    success = false;
    
    % Počiatočné odhady (vzdialené od seba, aby sme sa vyhli triviálnemu riešeniu)
    % R-fáza (bohatá na zložku 1), E-fáza (bohatá na zložku 2)
    xR = [0.98, 0.01, 0.01]; 
    xE = [0.01, 0.98, 0.01];
    
    diff_val = 1; 
    iter = 0;
    
    while diff_val > tol && iter < max_iter
        iter = iter + 1;
        
        gR = nrtl_gamma_calc(xR, delta_g, alpha, R, T);
        gE = nrtl_gamma_calc(xE, delta_g, alpha, R, T);
        
        K = gR ./ gE;
        xE_old = xE;
        
        % Successive substitution s tlmením (damping) pre stabilitu v kolóne
        xE_new = xE ./ sum(K .* xE);
        xE = 0.2 * xE_new + 0.8 * xE_old; 
        
        xR = K .* xE;
        xR = xR / sum(xR);
        
        diff_val = sum(abs(xE - xE_old));
    end

    % Kontrola stability (či sa fázy od seba líšia)
    if sum(abs(xR - xE)) > 1e-3
        success = true;
    end
end

function gamma = nrtl_gamma_calc(x, delta_g, alpha, R, T)
    n = length(x);
    tau = (delta_g .* 1000) ./ (R * T);
    G = exp(-alpha .* tau);
    S = x * G;
    
    term1 = (x * (tau .* G)) ./ S;
    term2 = zeros(1, n);
    for i = 1:n
        num_m = sum(x .* tau(:,i)' .* G(:,i)'); % Oprava indexácie pre konzistenciu
        % Pre ternárnu zmes v kolóne je bezpečnejšie použiť explicitnú sumu
        sum_j = 0;
        for j = 1:n
            num_m_j = sum(x .* tau(:,j)' .* G(:,j)');
            sum_j = sum_j + (x(j) * G(i,j) / S(j)) * (tau(i,j) - (num_m_j / S(j)));
        end
        term2(i) = sum_j;
    end
    gamma = exp(term1 + term2);
end