% Function to calculate the TE scattering cross section of a cylinder 
% 
% omega     = angular frequency [rad/s].
% a         = cylinder diameter [m].
% epsilon_d = permittivity of cylinder [1].
% epsilon_m = permittivity of surrounding medium [1].
% mu_m      = permeability of surrounding medium [1].
% mu_d      = permeability of cylinder [1].
%
% 
% HINT: make use of the besselj() and besselh() functions

function sigma_sc = scatteringCross(omega, a, epsilon_m, epsilon_d, mu_m, mu_d)
    c = 3e8; % Speed of light in a vacuum.
    k_m = (omega./c) * sqrt(epsilon_m * mu_m); % Wavenumber in medium
    k_d = (omega./c) * sqrt(epsilon_d * mu_d); % Wavenumber in cylinder
    
    % Radius of the cylinder
    R = a / 2;
    x = k_m * R; % Size parameter in medium
    y = k_d * R; % Size parameter in cylinder
    
    % Refractive index / Impedance ratios
    m = k_d / k_m; 
    mu_rel = mu_d / mu_m;
    
    % Determine number of terms for convergence
    N_max = ceil(x + 4 * x^(1/3) + 2);
    n = (0:N_max)';
    
    % 1. Calculate Bessel and Hankel functions
    % J_n(x), J_n(y), and H_n(1)(x)
    Jn_x = besselj(n, x);
    Jn_y = besselj(n, y);
    Hn_x = besselh(n, 1, x);
    
    % 2. Calculate derivatives using recurrence relations: 
    % Z_n'(z) = n/z * Z_n(z) - Z_{n+1}(z)  OR  0.5*(Z_{n-1}(z) - Z_{n+1}(z))
    % We use the identity: Z_n'(z) = Z_{n-1}(z) - (n/z)*Z_n(z)
    
    % For n=0, J_{-1}(z) = -J_1(z)
    Jn_prev_x = [-besselj(1, x); Jn_x(1:end-1)];
    Jn_prev_y = [-besselj(1, y); Jn_y(1:end-1)];
    Hn_prev_x = [-besselh(1, 1, x); Hn_x(1:end-1)];
    
    Jn_prime_x = Jn_prev_x - (n ./ x) .* Jn_x;
    Jn_prime_y = Jn_prev_y - (n ./ y) .* Jn_y;
    Hn_prime_x = Hn_prev_x - (n ./ x) .* Hn_x;
    
    % 3. Calculate TE coefficients (b_n)
    % Formula for TE: b_n = [ m*Jn'(y)*Jn(x) - mu_rel*Jn(y)*Jn'(x) ] / 
    %                       [ m*Jn'(y)*Hn(x) - mu_rel*Jn(y)*Hn'(x) ]
    % Note: m = k_d/k_m. For mu_rel = 1, this simplifies significantly.
    
    num = m * Jn_prime_y .* Jn_x - mu_rel * Jn_y .* Jn_prime_x;
    den = m * Jn_prime_y .* Hn_x - mu_rel * Jn_y .* Hn_prime_x;
    
    bn = num ./ den;
    
    % 4. Sum the series
    % sigma_sc = (4/k_m) * (|b0|^2 + 2 * sum_{n=1}^N |bn|^2)
    sum_terms = abs(bn(1))^2 + 2 * sum(abs(bn(2:end)).^2);
    sigma_sc = (4 / k_m) * sum_terms;
end
