% Function to solve for the odd modes of a symmetric slab waveguide 
% 
% ------------------
% >>>>> INPUTS <<<<<
% ------------------
% NOTE: Feel free to change what varibles this function accepts, you can
% add or remove anything depending on how you'd like to solve the problem.
%
% m = number mode to solve for [1].
% d = waveguide thickness in z direction [m].
% n1 = cladding refractive index
% n2 = core refractive index
% 
%
% ------------------
% <<<<< OUTPUT >>>>>
% ------------------
% omega = odd mode frequency [rad s-1].
% beta = odd mode wavevector (k_x) [m-1] 
% 

function [omega, beta] = findOddMode(m, d, n1, n2)
    % m: odd mode index (1, 3, 5...)
    % d: waveguide thickness [m]
    % n1: cladding index
    % n2: core index
    
    c = 3e8; % [m/s]
    
    % For the m-th odd mode, u resides in (m*pi/2, (m+1)*pi/2).
    % Sweep u across this range to generate the dispersion curve.
    u = linspace(m*pi/2 + 0.01, (m+1)*pi/2 - 0.01, 100);
    
    % 1. From the guidance condition for odd modes: v = -u * cot(u)
    v = -u ./ tan(u);
    
    % 2. Use the circle constraint to find the required V-parameter
    % V^2 = u^2 + v^2
    V = sqrt(u.^2 + v.^2);
    
    % 3. Relate V back to frequency (omega)
    % V = (omega * d / (2 * c)) * sqrt(n2^2 - n1^2)
    omega = (2 * V * c) ./ (d * sqrt(n2^2 - n1^2));
    
    % 4. Find the longitudinal wavevector (beta)
    % Using the separation of variables: k_core^2 = beta^2 + kz^2
    % (n2 * omega/c)^2 = beta^2 + (2u/d)^2
    k0 = omega ./ c;
    kz = 2 .* u ./ d;
    beta = sqrt((n2 .* k0).^2 - kz.^2);
    
end

