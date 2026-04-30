% Function to calculate the TE plane wave scattering field of a cylinder
% 
% r         = radius [m]
% phi       = rangle [rad]
% omega     = angular frequency [rad/s].
% a         = cylinder diameter [m].
% epsilon_d = permittivity of cylinder [1].
% epsilon_m = permittivity of surrounding medium [1].
% mu_m      = permeability of surrounding medium [1].
% mu_d      = permeability of cylinder [1].
%
% 
% HINT: make use of the besselj() and besselh() functions

function H = scatteringField(r, phi, omega, a, epsilon_m, epsilon_d, mu_m, mu_d)
    c = 3e8; % Speed of light in a vacuum.
    k_m = (omega./c) * sqrt(epsilon_m * mu_m); 
    k_d = (omega./c) * sqrt(epsilon_d * mu_d);
    
    R = a / 2;
    x = k_m * R; 
    y = k_d * R; 
    
    % Impedance and refractive index ratios
    m = k_d / k_m; 
    mu_rel = mu_d / mu_m;
    
    % Convergence criteria for the number of terms
    N_max = ceil(x + 4 * x^(1/3) + 2);
    n_vec = (-N_max:N_max)'; % Use symmetric indices for the field summation
    
    % 1. Calculate coefficients b_n for the cylinder
    % We only need unique n values for b_n since b_n = b_{-n}
    abs_n = abs(n_vec);
    
    Jn_x = besselj(abs_n, x);
    Jn_y = besselj(abs_n, y);
    Hn_x = besselh(abs_n, 1, x);
    
    % Derivatives using recurrence: Z_n' = Z_{n-1} - (n/z)*Z_n
    % Note: besselj(-1, z) = -besselj(1, z)
    Jn_prev_x = besselj(abs_n - 1, x);
    Jn_prev_y = besselj(abs_n - 1, y);
    Hn_prev_x = besselh(abs_n - 1, 1, x);
    
    Jn_prime_x = Jn_prev_x - (abs_n ./ x) .* Jn_x;
    Jn_prime_y = Jn_prev_y - (abs_n ./ y) .* Jn_y;
    Hn_prime_x = Hn_prev_x - (abs_n ./ x) .* Hn_x;
    
    % TE Scattering Coefficients b_n
    num = m * Jn_prime_y .* Jn_x - mu_rel * Jn_y .* Jn_prime_x;
    den = m * Jn_prime_y .* Hn_x - mu_rel * Jn_y .* Hn_prime_x;
    bn = num ./ den;
    
    % 2. Calculate the scattered field H_z at (r, phi)
    % H_z^s = sum [ bn * i^n * H_n(k_m*r) * exp(i*n*phi) ]
    H = 0;
    for i = 1:length(n_vec)
        n = n_vec(i);
        % Hankel function at the observation point r
        Hn_kr = besselh(n, 1, k_m * r);
        
        % Accumulate the field term
        term = bn(i) * (1i^n) * Hn_kr .* exp(1i * n * phi);
        H = H + term;
    end
end
