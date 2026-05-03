% Function to calculate the reflection coefficients of a planar
% air-dielectric interface.
% 
% omega = angular frequency.
% epsilon = dielectric function.
% theta_i = angle of incidence.
% r_p = p-polarized reflection coefficient.
% r_s = s-polarized reflection coefficient.
% 
function [r_p, r_s] = reflCoeff(omega, epsilon, theta_i)
    c = 3e8; % Speed of light in a vacuum.
    k_0 = omega/c;
    k_x = k_0*sin(theta_i*pi/180); % x-component of the incident wavevector.
    k_zi = sqrt(k_0^2 - k_x^2); % z-component of the incident wavevector.
    k_zt = sqrt(epsilon*k_0^2 - k_x^2); % z-component of the transmitted wavevector.
    
    % === TYPE YOUR CODE BELOW === %
    r_p = (k_zt - epsilon*k_zi)/(k_zt + epsilon*k_zi);
    r_s = (k_zi - k_zt)/(k_zi + k_zt);
end