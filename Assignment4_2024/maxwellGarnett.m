% Function to calculate the effective index of a porous material
% 
% ------------------
% >>>>> INPUTS <<<<<
% ------------------
% eps_m = metal permittivity [1] (vector or double)
% eps_a = ambient permittivity in pores [1] (vector or double)
% p = porosity, percent reported in decimal [1] (vector or double)
% ------------------
% <<<<< OUTPUT >>>>>
% ------------------
% eps_ff = effectric permittivity of the medium [1] (vector)
% 

function [eps_eff] = maxwellGarnett(eps_m,eps_a,p)
    % eps_m = metal permittivity
    % eps_a = ambient permittivity in pores
    % p = porosity (decimal, e.g., 0.1 for 10%)
    
    num = 2*eps_m + eps_a + 2*p.*(eps_a - eps_m);
    den = 2*eps_m + eps_a - p.*(eps_a - eps_m);
    
    eps_eff = eps_m .* (num ./ den);
end
        