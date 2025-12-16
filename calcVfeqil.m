function [alphaPF,betaPF] = calcVfeqil(T_farenheit, curve_name)

if ~ismember(curve_name, {'Castro', 'T27Aero', 'T27DT'})
    error('Not an implemented Beta Approach Curve')
end

switch curve_name
    case 'Castro'
        T_celsius = (T_farenheit - 32) * 5 / 9;
        % Castro Seraphin 1966 Beta Approach Curve
        betaPF = (7.5+92.5*exp(-0.0085*(980-T_celsius)))/100;
end

% Clamp the upper and lower bounds of the phase fraction...also, VPSC doesn't like it if
% you claim it's two-phase but have one phase set to a fraction of 0, so set it to close to 0
% instead.
betaPF(betaPF >= 1) = 0.999;
betaPF(betaPF <= 0) = 0.001;

alphaPF = 1-betaPF;

