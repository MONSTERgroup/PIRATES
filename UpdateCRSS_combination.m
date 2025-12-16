%% Function to update CRSS values

function [CRSS_basal, CRSS_prism, CRSS_pyra, CRSS_pyrca, CRSS_beta] = UpdateCRSS_combination(Tk)

%Li, Mason,  Bieler, Boehlert, Crimp, Acta Mat. 61(20)
%Table 1 (Reference 15)
basalp1 = [25+273.15 400];
prismp1 = [25+273.15 380];
pyramp1 = [25+273.15 640];

%Miller, Semiatin, Szczepanski, Pilchak, Met. Mat Trans A, 49(8), 2018  
% Table 5 (Reference 52)
% Beta ration ~= 1.2*Prism for Ti-6242 at room temperature, and that makes
% reasonable sense for Ti-64 as well.
betap1 = [25+273.15 1.2*prismp1(2)];

%From: Semiatin et al. MMT:A 33A 2002 p2719
%Lower bound for that equation's validity is 825C
basalp2 = [825+273.15 221.2];
prismp2 = [825+273.15 148.12];
pyramp2 = [825+273.15 253.15];
%Lower bound for beta equation's validity is 700C
betap2 = [700+273.15 109.9];

% Write current slip system files (temp. compensation)
R = 8.3143;
T0 = 1005.4;%K:  Temp where parameter ratios were determined. ==1350F

% Calculate temperature multipliers for the strengths of the alpha and
% beta phases. This is based on the formulation used in Semiatin et al
% MMT:A 33A 2002 p2719
Tmult_a = (exp((273000./(R.*Tk))-(273000/(R*T0)))).^(1/4.6);
Tmult_b = (exp((160000./(R.*Tk))-(160000/(R*T0)))).^(1/4.2);

%Miller, Semiatin, Szczepanski, Pilchak, Met. Mat Trans A, 49(8), 2018
RSSa0 = [269.8 402.9 461.1]; %Initial RSS values at 732C
RSSb0 = 94.5;

%Piecewise equation for alpha CRSS
if Tk > 825+273.15 + 1000000
    RSSa = RSSa0 * Tmult_a;
    
else
    RSSa(1) = linearize(prismp1, prismp2, Tk);
    RSSa(2) = linearize(basalp1, basalp2, Tk);
    RSSa(3) = linearize(pyramp1, pyramp2, Tk);
end

%Piecewise equation for beta CRSS
if Tk > 700+273.15
    RSSb = RSSb0 * Tmult_b;
else
    RSSb = linearize(betap1, betap2, Tk);
end

disp([RSSa RSSb]);

% Edit the *.SX files

CRSS_basal = RSSa(2);
CRSS_prism = RSSa(1);
CRSS_pyrca = RSSa(3);
CRSS_pyra = CRSS_prism;
CRSS_beta = RSSb;

end

function [y] = linearize(p1, p2, x)
y = (p2(2)-p1(2))/(p2(1)-p1(1))*(x-p1(1))+p1(2);
end