% Initialize VPSC for DEFORM-VPSC run

function [run] = initializeVPSCrun(path, csa, csb, vpsc_name, vpsc_path)

% This is where you can change the names of the input files for vpsc that you put in
% \Initial Inputs\

parameters = vpscParameters; 
    parameters.fromfile(fullfile(path, 'vpsc8.in'));
    parameters.iSave = 1;
    parameters.interactionType = [3; 10];
betaSX = vpscSingleCrystal; 
    betaSX.fromfile(fullfile(path, 'beta.sx'));
alphaSX = vpscSingleCrystal;    
    alphaSX.fromfile(fullfile(path, 'alpha.sx'));
betaTex = vpscTexture(csb, specimenSymmetry('-1')); 
    betaTex.fromfile(fullfile(path, 'beta.tex'));
alphaTex = vpscTexture(csa, specimenSymmetry('-1')); 
    alphaTex.fromfile(fullfile(path, 'alpha.tex'));

% Do a fake vpsc run to get an initial POSTMORT.OUT 
rolling_process = vpscDeformation; rolling_process.fromfile(fullfile(path, 'tiny_increment.in'));

run = vpscRun(parameters);
run.single_crystal{1} = betaSX;
run.single_crystal{2} = alphaSX;
run.texture_in{1} = betaTex;
run.texture_in{2} = alphaTex;

run.processes{1} = rolling_process;
run.vpsc_executable_name = vpsc_name;
run.magic_vpsc_box_path = vpsc_path;

run.write_inputs_to_magic_vpsc_box;
run.call_VPSC_executable;
