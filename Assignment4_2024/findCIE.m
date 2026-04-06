function [x, y] = findCIE(omega, R)
    c = 2.9979e8; 
    
    % Force omega and R to be column vectors to ensure compatible sizes
    omega = omega(:);
    R = R(:);
    
    % Convert angular frequency to wavelength in nanometers [cite: 523, 530]
    lambda = (2 * pi * c ./ omega) * 1e9; 

    % Selector function S(x, y, z) defined in Equation 4 [cite: 321, 322]
    S = @(val, gamma, delta) (val < 0) .* gamma + (val >= 0) .* delta;

    % Multi-lobe Gaussian Fits for the 1931 standard observer [cite: 315, 316]
    % Use coefficients from Table 1 [cite: 323, 324]
    
    % x_bar calculation
    t1_x = (lambda - 442.0) .* S(lambda - 442.0, 0.0624, 0.0374);
    t2_x = (lambda - 599.8) .* S(lambda - 599.8, 0.0264, 0.0323);
    t3_x = (lambda - 501.1) .* S(lambda - 501.1, 0.0490, 0.0382);
    x_fit = 0.362 * exp(-0.5 * t1_x.^2) + 1.056 * exp(-0.5 * t2_x.^2) - 0.065 * exp(-0.5 * t3_x.^2);

    % y_bar calculation
    t1_y = (lambda - 568.8) .* S(lambda - 568.8, 0.0213, 0.0247);
    t2_y = (lambda - 530.9) .* S(lambda - 530.9, 0.0613, 0.0322);
    y_fit = 0.821 * exp(-0.5 * t1_y.^2) + 0.286 * exp(-0.5 * t2_y.^2);

    % z_bar calculation
    t1_z = (lambda - 437.0) .* S(lambda - 437.0, 0.0845, 0.0278);
    t2_z = (lambda - 459.0) .* S(lambda - 459.0, 0.0385, 0.0725);
    z_fit = 1.217 * exp(-0.5 * t1_z.^2) + 0.681 * exp(-0.5 * t2_z.^2);

    % Integrate the products of reflectance and color matching functions [cite: 244, 245, 246]
    X = trapz(lambda, R .* x_fit);
    Y = trapz(lambda, R .* y_fit);
    Z = trapz(lambda, R .* z_fit);

    % Normalized chromaticity coordinates [cite: 251]
    total = X + Y + Z;
    x = X / total;
    y = Y / total;
end