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
    c = 2.99792e8;
    k_0 = omega/c;
    k_x = k_0*sin(theta);
    % Initialization: this should be Gamma_N.
    r = reflCoeff(omega, k_x, epsilon(N+1), epsilon(N+2), pol);
    lastGamma = r; % lastGamma always stores the previous value of Gamma.
    
    % Loop over the remaining layers.
    for i = (N+1):(-1):2
        % === TYPE YOUR CODE BELOW === %
        
        
        
        
    end
    
    % Reflectance.
    % === TYPE YOUR CODE BELOW === %
    
end