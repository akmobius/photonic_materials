% Function to calculate the CIE x,y chromaticity coordinates from a reflectance
% spectrum
% 
% ------------------
% >>>>> INPUTS <<<<<
% ------------------
% omega = angular frequency [rad/s].
% R = reflectance [1].
% 
% ------------------
% <<<<< OUTPUT >>>>>
% ------------------
% x,y = CIE chromaticity coordinates [1].
% 
function [x, y] = findCIE(omega, R)
    % Force inputs to be row vectors (1xN) regardless of original shape
    omega = omega(:).'; 
    R = R(:).';
    
    c = 3e8;
    
    % Calculate wavelength in nm and force to 1xN
    lambda_nm = (2 * pi * c ./ omega) * 1e9;
    lambda_nm = lambda_nm(:).'; 

    % --- ANALYTIC APPROXIMATIONS (CIE 1931) ---
    % We apply ( : ).' to every analytic term to guarantee row vectors
    
    % X-bar: forced row
    x_bar = (1.056 * exp(-0.5 * ((lambda_nm - 599.8) / 37.9).^2) + ...
             0.362 * exp(-0.5 * ((lambda_nm - 442.0) / 16.0).^2) - ...
             0.065 * exp(-0.5 * ((lambda_nm - 501.1) / 20.4).^2));
    x_bar = x_bar(:).';

    % Y-bar: forced row
    y_bar = (0.821 * exp(-0.5 * ((lambda_nm - 568.8) / 46.9).^2) + ...
             0.286 * exp(-0.5 * ((lambda_nm - 530.9) / 16.3).^2));
    y_bar = y_bar(:).';

    % Z-bar: forced row
    z_bar = (1.217 * exp(-0.5 * ((lambda_nm - 437.0) / 22.6).^2));
    z_bar = z_bar(:).';

    % --- INTEGRATION ---
    % Since R and x_bar are strictly forced to rows, .* will not fail
    X = trapz(lambda_nm, R .* x_bar);
    Y = trapz(lambda_nm, R .* y_bar);
    Z = trapz(lambda_nm, R .* z_bar);
    
    % Normalization to find x, y coordinates
    total_sum = X + Y + Z;
    
    % Handle potential edge cases where R might be zero
    if total_sum == 0
        x = 0.3333; 
        y = 0.3333;
    else
        x = X / total_sum;
        y = Y / total_sum;
    end
end