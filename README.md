# PIRATES <img width="50" alt="Cartoon image of a PIRATE." class="recess" src="https://github.com/user-attachments/assets/a79ceb6f-7382-442a-af2b-e43df6961a03">
Code for the Process-Scale Industry-Research Automated Texture Evolution Simulation (PIRATES).


## Cite As
B.A. Begley, M.E. Hurley, M.S. David & V.M. Miller. "A Multiscale Model Predicting Site-Specific Texture Evolution: Application to Two-Phase Titanium Alloys." *Integrating Materials and Manufacturing Innovation* (Accepted). [DOI](https://doi.org/https://doi.org/10.1007/s40192-025-00433-2).

## Requirements and Dependencies

**MATLAB** (R2019a or later; earlier version may work but are untested)
- Some dependencies (MTEX, VAMPYR) use the MATLAB **Optimization Toolbox** for additional features, but this toolbox is not necessary for PIRATES to be used as demonstrated in the manuscript.

[**MTEX**](https://mtex-toolbox.github.io/download) open source MATLAB toolbox for crystallographic texture analysis
- **Version 5.8.1** is a currently a hard version requirement. Version 5.10+ changed syntax that breaks compatibility with PIRATES, and I have not personally tested PIRATES with versions 5.8.2 or 5.9.X.

[**VAMPYR**](https://github.com/MONSTERgroup/VAMPYR/tree/vampyr8): **V**PSC **A**utomation in **M**TEX for **P**olycrystal plasticit**Y** **R**esearch
- Specifically, the linked vampyr8 branch which has syntax for VPSC 8, rather than the initial version for VPSC 7.

[**VPSC 8**](https://github.com/lanl/VPSC_code) The viscoplastic self-consistent polycrystal plasticity model
- A VPSC executable compiled for Windows 11 is included included with VAMPYR. Users with incompatible systems will have to compile VPSC from source code.

## Initial Setup and Startup

In order to use PIRATES use the `startup_mtex()` function in your MTEX 5.8.1 directory (or install MTEX permanently to your MATLAB path) then use the `startup_api()` function in your VAMPYR directory. Navigate to the PIRATES directory.

PIRATES has one hardcoded path declaration in the [Main3](Main3.m) script, which must direct PIRATES to the "magic_vpsc_box" directory in your VAMPYR installation.

```matlab
% Define the location of the VAMPYR magic_vpsc_box and executable name.
vpsc_name = 'vpsc8.exe';
vpsc_path = 'C:\path_to_vampyr_installation\VPSC\magic_vpsc_box';
```

## Setting up a run

The input files should be placed in [Initial Inputs](Initial%20Inputs). This directory should include:
- [alpha.sx](Initial%20Inputs/alpha.sx) and [beta.sx](Initial%20Inputs/beta.sx). These are templates for the VPSC single crystal parameter files with pre-defined slip systems for &alpha; and &beta; Ti.
- [alpha.tex](Initial%20Inputs/alpha.tex) and [beta.texx](Initial%20Inputs/beta.tex). VPSC texture input files.
- [vpsc8.in](Initial%20Inputs/vpsc8.in). A template file for the main VPSC8 input file.
- [tiny_increment.in](Initial%20Inputs/tiny_increment.in). A template file for a VPSC deformation step.
- DEFORM Point Tracking Files as the .CSV file type.

## Pre-Run Checklist
### Settings that can be changed in [Main3](Main3.m)
These parameters control the simulation output.

#### Set the name of the point track file:
```matlab
% Name of the DEFORM point tracking file (Which should be located in:
%   \Initial Inputs\)
pointTrackFile = 'DEF_PTR001_1600_2p.CSV';
```
#### Choose to apply a random recrystallization texture or an &alpha; BOR variant texture at specific points:
```matlab
% Apply a random recrystallization texture to both phases at Reset Points?
isRandomRX = false; % If this is true, set isAlphaRandomVariant to false.

% Apply a random recrystallization texture to the beta phase and a BOR texture with
% random Alpha variant selection to the alpha phase at reset points?
isAlphaRandomVariant = false; % If this is true, set isRandomRX to false.
```
#### Choose to which &beta; approach curve to use to calculate phase fractions:
- Currently, only the curve from Castro & Seraphin (1966) is implemented.
- More methods can be added to the `switch` statement in [calcVfeqil](calcVfeqil.m).
```matlab
% Choose from: 'Castro';
beta_approach_curve = 'Castro';
```
```matlab
if ~ismember(curve_name, {'Castro'})
    error('Not an implemented Beta Approach Curve')
end

switch curve_name
    case 'Castro'
        T_celsius = (T_farenheit - 32) * 5 / 9;
        % Castro Seraphin 1966 Beta Approach Curve
        betaPF = (7.5+92.5*exp(-0.0085*(980-T_celsius)))/100;
end
```
#### Choose which steps recrystallization textures should be applied to (if any):
- Set `reset_step = [];` if it will be unused.
```matlab
% =========================================================================
% --- RESET STEPS ---
% This is a new feature!
% Basically, these are steps where we are going to choose to manually override the
% texture in the VPSC input files, in order to apply (for instance) a recrystallization
% texture. This can't be done at every step, because this also zeros out the residual
% stresses in VPSC, which are necessary for accurate plastic rotation prediction.
% My recommendation is to set the steps to be at the end of a long isothermal hold.
% =========================================================================

reset_steps = [129 265 456]; %1400 4p
%reset_steps = [129 265 456]; %1600 4p
%reset_steps = [];
```
#### Choose the desired `UpdateCRSS()` function:
- `UpdateCRSS_powerlaw()` is the method used for the vast majority of the manuscript.
- `UpdateCRSS_combined()` was used to demonstrate non-fixed CRSS ratios.
- `UpdateCRSS_linear()` is a slightly simplified version of `UpdateCRSS_combined()`.
```matlab
[CRSS_basal, CRSS_prism, CRSS_pyra, CRSS_pyrca, CRSS_beta] = UpdateCRSS_powerlaw(Tk);
run.single_crystal{1,2}.Modes{1,1}.voceParams(1,1) = CRSS_prism;
run.single_crystal{1,2}.Modes{2,1}.voceParams(1,1) = CRSS_basal;
run.single_crystal{1,2}.Modes{3,1}.voceParams(1,1) = CRSS_pyra;
run.single_crystal{1,2}.Modes{4,1}.voceParams(1,1) = CRSS_pyrca;
run.single_crystal{1,1}.Modes{1,1}.voceParams(1,1) = CRSS_beta;
run.single_crystal{1,1}.Modes{2,1}.voceParams(1,1) = CRSS_beta;
run.single_crystal{1,1}.Modes{3,1}.voceParams(1,1) = CRSS_beta;
```
#### Name the output directory:
```matlab
output_directory = 'outputs\1600_2p';
mkdir(fullfile(MainDir,output_directory));
```

### Settings that can be changed in [Analysis3](Analysis3.m):
These parameters control the plotting of figures from the simulation data.

#### Set the name of the point track file and output directory:
- These should be the same as in Main3.
```matlab
pointTrackFile = 'DEF_PTR001_1600_2p.CSV';
folder = 'outputs\1600_2p';
```
#### Set the list of tracked points to analyze:
```matlab
% List of points you want to analyze. Does not need to be the full list of points.
wantedPt = [1 2 4];
```
#### Give the points descriptive names:
```matlab
% Names of points for labeling figures. Needs to be as long as the number of points in
% the point tracking file.
pos_name = {'ND_Offset', 'Center', 'TD_Offset', 'RD_Offset'};
```
#### Make sure the crystal symmetry parameters are the same as in Main3:
```matlab
CSa = crystalSymmetry('622', [2.95 2.95 4.68], 'X||a', 'Y||b*', 'Z||c', 'color', 'light blue');
CSb = crystalSymmetry('432', [3.24 3.24 3.24], 'mineral', 'Titanium - Beta', 'color', 'light green');
```

## Run [Main3](Main3.m)
By default, [Analysis3](Analysis3.m) is called at the end of [Main3](Main3.m), so plotting will run automatically.



