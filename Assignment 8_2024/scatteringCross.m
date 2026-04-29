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
    k_m =  (omega./c)*sqrt(epsilon_m*mu_m); % m-1
    k_d =  (omega./c)*sqrt(epsilon_d*mu_d); % m-1
    
    % === TYPE YOUR CODE BELOW === % 


end