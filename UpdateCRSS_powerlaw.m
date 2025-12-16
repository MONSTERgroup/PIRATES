%% Function to update CRSS and TRSS values for Ti-64
function [CRSS_basal, CRSS_prism, CRSS_pyra, CRSS_pyrca, CRSS_beta] = UpdateCRSS_powerlaw(Tk)

% Effectively, find an equation for one CRSS, and assign the rest as a function of
% that.
CRSS_prism = 56.33 + 420.9 .* exp(-0.004453 .* Tk);
CRSS_basal = 1.5 .* CRSS_prism;
CRSS_pyra  = 1 .* CRSS_prism;
CRSS_pyrca  = 3 .* CRSS_prism;
CRSS_beta = 0.3 .* CRSS_prism;

end