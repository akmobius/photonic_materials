% Function to calculate the reflection coefficients of a planar
% interface between isotropic materials.
% 
% omega = angular frequency [rad/s].
% k_x = x-component of the wavevector [1/m].
% epsilon_1 = dielectric function of material 1, the incident medium,
%             evaluated at omega [1].
% epsilon_2 = dielectric function of material 2, the transmitting medium,
%             evaluated at omega [1].
% pol = polarization (p or s).
% 
function r = reflCoeff(omega, k_x, epsilon_1, epsilon_2, pol)
    c = 3e8; % Speed of light in a vacuum.
    k_0 = omega/c;
    k_z1 = sqrt(epsilon_1*k_0^2 - k_x^2); % z-component of the incident wavevector.
    k_z2 = sqrt(epsilon_2*k_0^2 - k_x^2); % z-component of the transmitted wavevector.
    
    % Calculate the reflection coefficient.
    if pol == 'p'
        r = (epsilon_1*k_z2 - epsilon_2*k_z1)/(epsilon_1*k_z2 + epsilon_2*k_z1);
    elseif pol == 's'
        r = (k_z2 - k_z1)/(k_z2 + k_z1);
    end
end