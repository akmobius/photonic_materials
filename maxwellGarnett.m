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
        eps_eff = eps_m .* ( (2.*(1-p).*eps_m + (1+2.*p).*eps_a) ./ ((2+p).*eps_m + (1-p).*eps_a) );
        
end