%% Function to update CRSS and TRSS values for Ti-64
function [CRSS_basal, CRSS_prism, CRSS_pyra, CRSS_pyrca, CRSS_beta] = UpdateCRSS_linear(Tk)

%% Define low and high temperature data points (MPa)

% [Temperature(K), Stress(MPa)]
bounds = struct(...
    'basal',  [298.15, 400; 1098.15, 221.2], ...
    'prism',  [298.15, 380; 1098.15, 148.12], ...
    'pyra',   [298.15, 380+10; 1098.15, 148.12+10], ...% Assume pyr<a> ~= prism (but with a small offset).
    'pyrca',  [298.15, 640; 1098.15, 253.15], ...
    'beta',   [298.15, 456; 973.15, 109.9]);  % Beta ~1.2*prism at RT

% Alpha phase
CRSS_basal = interpolate(bounds.basal, Tk);
CRSS_prism = interpolate(bounds.prism, Tk);
CRSS_pyra  = interpolate(bounds.pyra, Tk);
CRSS_pyrca = interpolate(bounds.pyrca, Tk);

% Beta phase
CRSS_beta = interpolate(bounds.beta, Tk);
% Fast way of making it so that the value can't go below zero within a range
CRSS_beta(Tk > 1150) = CRSS_beta(Tk==1150);

end

%% Helper function for linear interpolation
function y = interpolate(points, x)
% points: 2x2 matrix [T1, stress1; T2, stress2]
x1 = points(1,1);
x2 = points(2,1);
y1 = points(1,2);
y2 = points(2,2);

y = y1 + (y2 - y1) ./ (x2 - x1) .* (x - x1);
end