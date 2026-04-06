% Function to calculate the reflectance of any multilayer structure.
% 
% ------------------
% >>>>> INPUTS <<<<<
% ------------------
% omega = angular frequency [rad/s].
% theta = polar angle of incidence [rad].
% N = number of layers [1].
% epsilon = (N+2)x1 array containing the dielectric function of each layer
%           evaluated at omega. The first element should correspond to the
%           0th layer; the last element should correspond to the (N+1)th
%           layer.
% d = Nx1 array containing the thickness of each layer, except the 0th and
%     (N+1)th layer, which are semi-infinite.
% pol = polarization (p or s).
% 
% ------------------
% <<<<< OUTPUT >>>>>
% ------------------
% R = reflectance [1].
% 
function R = multilayer(omega, theta, N, epsilon, d, pol)
    % Setup.
    c = 2.99792e8; % 
    k_0 = omega/c; % 
    k_x = k_0*sin(theta); % 

    % Initialization: Calculate Gamma at the very last interface.
    % In a system with N layers, the last interface is between 
    % layer N (index N+1) and the substrate (index N+2). 
    lastGamma = reflCoeff(omega, k_x, epsilon(N+1), epsilon(N+2), pol); 
    
    % Loop over the remaining layers from the substrate back to the incident medium.
    % We iterate from the last physical layer (N) up to the first physical layer (1). 
    for i = N:(-1):1
        % 1. Calculate the vertical wavevector in the current layer i.
        % Layer i corresponds to epsilon(i+1) because epsilon(1) is the incident medium.
        k_zi = sqrt(epsilon(i+1)*k_0^2 - k_x^2);
        
        % 2. Calculate the phase factor (propagation) across layer i of thickness d(i).
        phase = exp(2j * k_zi * d(i)); 
        
        % 3. Calculate the local Fresnel reflection coefficient at the interface 
        % between layer i-1 (epsilon(i)) and layer i (epsilon(i+1)).
        r_local = reflCoeff(omega, k_x, epsilon(i), epsilon(i+1), pol); 
        
        % 4. Apply the recursive formula for the generalized reflection coefficient.
        lastGamma = (r_local + lastGamma * phase) / (1 + r_local * lastGamma * phase); 
    end
    
    % The final reflectance is the square of the magnitude of the 
    % reflection coefficient at the first interface.
    R = abs(lastGamma)^2;