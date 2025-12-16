% =========================================================================
% --- PIRATES ---
% Process-Scale Industry-Research Automated Texture Evolution Simulation
% As demonstrated in: "A Multiscale Model Predicting Site-Specific Texture
% Evolution: Application to Two-Phase Titanium Alloys."
% =========================================================================

% start timing the run (left over from debug, but sill useful)
tic

% =========================================================================
% --- RUN SETUP ---
% =========================================================================

% This should be the directory that Main3.m is in.
MainDir = pwd;
if ~exist(fullfile(pwd, 'Initial Inputs'), 'dir')
    error('Cannot find initial inputs folder; confirm that the script is being called from the correct working directory.');
end

% Define the location of the VAMPYR magic_vpsc_box and executable name.
vpsc_name = 'vpsc8.exe';
vpsc_path = 'E:\Final_PIRATES\For_GitHub\vampyr\VPSC\magic_vpsc_box';

% Name of the DEFORM point tracking file (Which should be located in:
%   \Initial Inputs\)
pointTrackFile = 'DEF_PTR001_1600_2p.CSV';

% This parameter is left over from testing the framework. Replace T_fahrenheight with
% this parameter if you want to assume a fixed temperature instead of using the
% temperature from the DEFORM simulation
nominal_T_fahrenheit = 1400;

% Apply a random recrystallization texture to both phases at Reset Points?
isRandomRX = false; % If this is true, set isAlphaRandomVariant to false.

% Apply a random recrystallization texture to the beta phase and a BOR texture with
% random Alpha variant selection to the alpha phase at reset points?
isAlphaRandomVariant = false; % If this is true, set isRandomRX to false.

% Choose from: 'Castro';
beta_approach_curve = 'Castro';

% Print extra information to the Command Window?
verbose = true;

% =========================================================================
% --- READ DATA FROM POINT TRACK FILE ---
% =========================================================================

% Read and segment the point track file
[ptList,nIter] = DeterminePointsIter(fullfile(MainDir,'Initial Inputs',pointTrackFile));
[vel_grad,strain_inc,temp, time, pos,step, strain_rate] = readDEFPTR3d(fullfile(MainDir,'Initial Inputs',pointTrackFile),length(ptList),nIter);

% =========================================================================
% --- MAYBE HELPFUL FIGURES ---
% These figures plot the position of points in the point track file relative to
% each other for (1) all steps and (2) just the first step.
% =========================================================================

figure;
scatter3(pos(:,1,1),pos(:,1,2),pos(:,1,3))
hold on
scatter3(pos(:,2,1),pos(:,2,2),pos(:,2,3))
scatter3(pos(:,3,1),pos(:,3,2),pos(:,3,3))
scatter3(pos(:,4,1),pos(:,4,2),pos(:,4,3))
hold off
xlabel('X')
ylabel('Y')
zlabel('Z')
legend({'1','2','3','4'})
PrettyPlotsSingle;

figure;
scatter3(pos(1,1,1),pos(1,1,2),pos(1,1,3),'filled')
hold on
scatter3(pos(1,2,1),pos(1,2,2),pos(1,2,3),'filled')
scatter3(pos(1,3,1),pos(1,3,2),pos(1,3,3),'filled')
scatter3(pos(1,4,1),pos(1,4,2),pos(1,4,3),'filled')
hold off
xlabel('X')
ylabel('Y')
zlabel('Z')
legend({'1','2','3','4'})
PrettyPlotsSingle;

%%

% =========================================================================
% --- A BIT MORE RUN PREPARATION ---
% =========================================================================

% Segment the point track file (into deformation, hold, etc. steps)
[~, seg_mat] = SegmentPointTrack(fullfile(MainDir,'Initial Inputs',pointTrackFile),pos);

% Parameters for calculating ODFs
% I've changed these to match the parameters from the .ctf file I received!
CSa = crystalSymmetry('622', [2.95 2.95 4.68], 'X||a', 'Y||b*', 'Z||c', 'color', 'light blue');
CSb = crystalSymmetry('432', [3.24 3.24 3.24], 'mineral', 'Titanium - Beta', 'color', 'light green');
SS = specimenSymmetry('-1');
psi = deLaValleePoussinKernel('halfwidth', 8*degree);

% Track how many "deformation" steps get removed because the deformation increment is
% too small for VPSC to handle (this mostly just catches floating point errors that
% should be zeros)
skipped = zeros(1,4);


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

%%

% =========================================================================
% --- LOOP OVER POINTS in the POINT TRACKING FILE ---
% =========================================================================

% if you want to do all points, use:
% for pt = ptList'

% manually choose the points you want to simulate out of the file.
for pt = [2]

    % Set up empty output variables
    %   Phase fraction [beta alpha]
    phaseFraction = [ones(nIter,1) zeros(nIter,1)];

    % ODFs
    outputODF = [repmat(uniformODF(CSa,SS),nIter,1) repmat(uniformODF(CSb,SS),nIter,1)];
    cumulativeRot = repmat(rotation.byAxisAngle(xvector, 0*degree),nIter,1);

    activities = struct();
    activities.prism = zeros(nIter, 1);
    activities.basal = zeros(nIter, 1);
    activities.pyra  = zeros(nIter, 1);
    activities.pyrca = zeros(nIter, 1);

    % get just the deformation quantities for the current point
    [vel_pt,strain_inc_pt,temp_pt,time_pt,step_pt, sr_pt] = subsetDefPTR(pt,time,temp,strain_inc,vel_grad,step,strain_rate);

    % Set up initial VPSC run
    run = initializeVPSCrun('Initial Inputs', CSa, CSb, vpsc_name, vpsc_path);
    old_run = run;

    % Make future runs look for a POSTMORT file.
    old_run.parameters.iRecover = 1;
    old_run.parameters.iSave = 11;

    % Set the phase fractions and ODFs for the first step.
    phaseFraction(1,:) = old_run.parameters.phaseFrac;
    outputODF(1,1) = calcDensity(old_run.texture_out{1}.orientations,'weights',old_run.texture_out{1}.weights, 'kernel', psi);
    outputODF(1,2) = calcDensity(old_run.texture_out{2}.orientations,'weights',old_run.texture_out{2}.weights, 'kernel', psi);

    activities.prism(1) = 0;
    activities.basal(1) = 0;
    activities.pyra(1)  = 0;
    activities.pyrca(1) = 0;

    % =========================================================================
    % --- LOOP OVER DEFORM STEPS for each POINT ---
    % =========================================================================

    for it = 2:nIter

        run = vpscRun(old_run,'follow');

        [curr_vel,curr_inc,delta_T, curr_sr] = getIncrVals(it,vel_pt,strain_inc_pt,temp_pt,sr_pt);
        T_fahrenheit = temp_pt(it);

        % =========================================================================
        % --- LOGIC for RESET POINTS ---
        % =========================================================================

        if ismember(it, reset_steps)
            % Use a different method of setting up the VPSC run.
            run = vpscRun(old_run,'follow_no_postmort');

            % Update phase frac
            [alphaPF,betaPF] = calcVfeqil(T_fahrenheit, beta_approach_curve);
            run.parameters.phaseFrac = [betaPF alphaPF];
            phaseFraction(it,:) = run.parameters.phaseFrac;

            run.parameters.iRecover = 0;

            % do a fake deformation process to generate a POSTMORT
            run.parameters.iSave = 1;
            run.parameters.nProcess = 1;
            run.parameters.processType = 0;
            rolling_process = vpscDeformation; rolling_process.fromfile(fullfile('Initial Inputs', 'tiny_increment.in'));
            run.processes{1} = rolling_process;

            if isRandomRX % Logic for random recrystallization texture

                temp_b = calcDensity(run.texture_in{1}.orientations,'weights',run.texture_in{1}.weights, 'kernel', psi);
                temp_a = calcDensity(run.texture_in{2}.orientations,'weights',run.texture_in{2}.weights, 'kernel', psi);

                % Make a weighted average of the VPSC input texture and a uniform
                % texture. Use these parameters to strengthen or weaken each component.
                temp_b = 0.8*temp_b + 0.2*uniformODF(CSb,SS);
                temp_a = 0.8*temp_a + 0.2*uniformODF(CSa,SS);

                temp_b = temp_b.discreteSample(run.texture_in{1}.ngrain);
                temp_a = temp_a.discreteSample(run.texture_in{2}.ngrain);

                run.texture_in{1}.orientations = temp_b;
                run.texture_in{2}.orientations = temp_a;

            elseif isAlphaRandomVariant % (Naive) Logic for BOR texture for Alpha phase
                % Basically just uses MTEX parent_child relationship to generate
                % possible alpha orientation distributions.
                % could be updated to include variant selection rules, somehow.

                beta2alpha = orientation.Burgers(CSb,CSa);

                oriParent = run.texture_in{1}.orientations;
                oriChild = variants(beta2alpha, oriParent);

                temp_b = calcDensity(run.texture_in{1}.orientations,'weights',run.texture_in{1}.weights, 'kernel', psi);
                temp_a = calcDensity(run.texture_in{2}.orientations,'weights',run.texture_in{2}.weights, 'kernel', psi);
                temp_a_child = calcDensity(oriChild, 'kernel', psi);

                temp_b = 0.8*temp_b + 0.2*uniformODF(CSb,SS);
                temp_a = 0.8*temp_a + 0.2*temp_a_child;

                temp_b = temp_b.discreteSample(run.texture_in{1}.ngrain);
                temp_a = temp_a.discreteSample(run.texture_in{2}.ngrain);

                run.texture_in{1}.orientations = temp_b;
                run.texture_in{2}.orientations = temp_a;

            else
                % if both are set to false, do nothing at the reset points. That said,
                % you probably want to empty the reset_points array if that's how it's
                % set up, then you'll never get here in the first place.
            end


            % Run VPSC
            run.vpsc_executable_name = vpsc_name;
            run.magic_vpsc_box_path = vpsc_path;
            run.write_inputs_to_magic_vpsc_box;
            run.call_VPSC_executable;
            run.cleanup;
            old_run = run;

            % Set it back to normal iSave values.
            old_run.parameters.iSave = 11;
            old_run.parameters.nProcess = 1;
            old_run.parameters.processType = 0;
            %.parameters.processDetail(2,:) = [];
        end

        % =========================================================================
        % --- LOGIC for the 4 types of steps ---
        % =========================================================================

        switch seg_mat{it}
            case 'nothing'
                % i.e., it's a remesh or something else that doesn't affect the
                % simulation in a physical way.
                if verbose; fprintf('Step: %i, doing: %s\n', it, seg_mat{it}); end

                % Copy phase frac and ODF from previous step
                phaseFraction(it,:) = old_run.parameters.phaseFrac;
                outputODF(it,1) = calcDensity(old_run.texture_out{1}.orientations,'weights',old_run.texture_out{1}.weights, 'kernel', psi);
                outputODF(it,2) = calcDensity(old_run.texture_out{2}.orientations,'weights',old_run.texture_out{2}.weights, 'kernel', psi);
                cumulativeRot(it) = cumulativeRot(it-1);

                activities.prism(it) = 0;
                activities.basal(it) = 0;
                activities.pyra(it)  = 0;
                activities.pyrca(it) = 0;

            case 'hold'
                % i.e., an isothermal hold; the temperature changed, but there's no
                % deformation at the point.
                if verbose; fprintf('Step: %i, doing: %s\n', it, seg_mat{it}); end

                % Update phase frac, copy ODF from prev step
                [alphaPF,betaPF] = calcVfeqil(T_fahrenheit, beta_approach_curve);
                old_run.parameters.phaseFrac = [betaPF alphaPF];
                phaseFraction(it,:) = old_run.parameters.phaseFrac;

                outputODF(it,1) = calcDensity(old_run.texture_out{1}.orientations,'weights',old_run.texture_out{1}.weights, 'kernel', psi);
                outputODF(it,2) = calcDensity(old_run.texture_out{2}.orientations,'weights',old_run.texture_out{2}.weights, 'kernel', psi);
                cumulativeRot(it) = cumulativeRot(it-1);

                activities.prism(it) = 0;
                activities.basal(it) = 0;
                activities.pyra(it)  = 0;
                activities.pyrca(it) = 0;

            case 'deformation'
                % i.e., VPSC actually needs to run.
                if verbose; fprintf('Step: %i, doing: %s\n', it, seg_mat{it}); end

                run = vpscRun(old_run,'follow');

                % Update phase frac
                [alphaPF,betaPF] = calcVfeqil(T_fahrenheit, beta_approach_curve);
                run.parameters.phaseFrac = [betaPF alphaPF];
                phaseFraction(it,:) = run.parameters.phaseFrac;

                % =========================================================
                % --- NOTE ---
                % I've changed it so it only uses variants of UpdateCRSS()
                % instead of the original method. Change the function name
                % here (UpdateCRSS_linear, UpdateCRSS_powerlaw,
                % UdateCRSS_yourcustomfile) to change how CRSSs are
                % calculated.


                Tk = (T_fahrenheit + 459.67) * (5/9);
                [CRSS_basal, CRSS_prism, CRSS_pyra, CRSS_pyrca, CRSS_beta] = ...
                    UpdateCRSS_powerlaw(Tk);
                run.single_crystal{1,2}.Modes{1,1}.voceParams(1,1) = CRSS_prism;
                run.single_crystal{1,2}.Modes{2,1}.voceParams(1,1) = CRSS_basal;
                run.single_crystal{1,2}.Modes{3,1}.voceParams(1,1) = CRSS_pyra;
                run.single_crystal{1,2}.Modes{4,1}.voceParams(1,1) = CRSS_pyrca;
                run.single_crystal{1,1}.Modes{1,1}.voceParams(1,1) = CRSS_beta;
                run.single_crystal{1,1}.Modes{2,1}.voceParams(1,1) = CRSS_beta;
                run.single_crystal{1,1}.Modes{3,1}.voceParams(1,1) = CRSS_beta;

                % =========================================================

                % Update deformation process
                def_process = vpscDeformation; def_process.fromfile(fullfile('Initial Inputs', 'tiny_increment.in'));

                % Uncomment this line to use the velocity gradient
                %L =  [curr_vel(1:3);curr_vel(4:6);curr_vel(7:9)];

                % This line uses the symmetric strain rate tensor instead.
                L =  [curr_sr(1:3);curr_sr(4:6);curr_sr(7:9)];

                % These values seem to work well.
                def_process.nsteps = 10;
                def_process.ictrl = 3;
                def_process.increment = abs(curr_inc(3))/def_process.nsteps;

                % DEFORM doesn't always give you tensors where the trace of the matrix
                % is zero, but VPSC requires that. This distributes the error over the
                % entire diagonal in order to force the trace to be zero. It's ugly, but
                % it works.
                if trace(L) ~= 0
                    shift_amount = trace(L)/3;
                    L(1) = L(1)-shift_amount;
                    L(5) = L(5)-shift_amount;
                    L(9) = L(9)-shift_amount;
                end

                skip_flag = false;
                if L(9) < 1e-8
                    [M,I] = max(abs([L(1) L(5)]),[],'all');
                    def_process.ictrl = I;
                    if M < 1e-8
                        skip_flag = true;
                    end
                end

                def_process.velocity_gradient = L;

                if def_process.increment < 1e-8 || skip_flag || not(any([L(1) L(5) L(9)],'all'))
                    %Fake a hold if the increment is too small
                    skipped(pt) = skipped(pt) + 1;
                    old_run.parameters.phaseFrac = [betaPF alphaPF];
                    phaseFraction(it,:) = old_run.parameters.phaseFrac;
                    outputODF(it,1) = calcDensity(old_run.texture_out{1}.orientations,'weights',old_run.texture_out{1}.weights, 'kernel', psi);
                    outputODF(it,2) = calcDensity(old_run.texture_out{2}.orientations,'weights',old_run.texture_out{2}.weights, 'kernel', psi);
                    cumulativeRot(it) = cumulativeRot(it-1);

                    activities.prism(it) = 0;
                    activities.basal(it) = 0;
                    activities.pyra(it)  = 0;
                    activities.pyrca(it) = 0;
                else
                    run.parameters.processType = 0;
                    run.processes{1} = def_process;

                    % Run VPSC
                    run.vpsc_executable_name = vpsc_name;
                    run.magic_vpsc_box_path = vpsc_path;
                    run.write_inputs_to_magic_vpsc_box;
                    run.call_VPSC_executable;
                    run.cleanup;
                    old_run = run;

                    % Get ODF from VPSC run
                    outputODF(it,1) = calcDensity(old_run.texture_out{1}.orientations,'weights',old_run.texture_out{1}.weights, 'kernel', psi);
                    outputODF(it,2) = calcDensity(old_run.texture_out{2}.orientations,'weights',old_run.texture_out{2}.weights, 'kernel', psi);
                    cumulativeRot(it) = cumulativeRot(it-1);

                    activities.prism(it) = old_run.slip_activity{2}.activities(end,1);
                    activities.basal(it) = old_run.slip_activity{2}.activities(end,2);
                    activities.pyra(it)  = old_run.slip_activity{2}.activities(end,3);
                    activities.pyrca(it) = old_run.slip_activity{2}.activities(end,4);
                end


            case 'rotation'

                run = vpscRun(old_run,'follow');
                % Update phase frac
                [alphaPF,betaPF] = calcVfeqil(T_fahrenheit, beta_approach_curve);
                run.parameters.phaseFrac = [betaPF alphaPF];
                phaseFraction(it,:) = run.parameters.phaseFrac;

                % Update rotation process
                rot_process = vpscRotation;
                rot_mat = calcRot(it,squeeze(pos(:,:,:)));
                rot_process.rot = rotation.byMatrix(rot_mat);
                run.processes{1} = rot_process;

                % do a fake deformation process to generate a POSTMORT
                run.parameters.iSave = 1;
                run.parameters.nProcess = 2;
                run.parameters.processType = [4; 0];
                run.parameters.processDetail{2,1} = 'dummy.def';
                rolling_process = vpscDeformation; rolling_process.fromfile(fullfile('Initial Inputs', 'tiny_increment.in'));
                run.processes{2} = rolling_process;

                % Run VPSC
                run.vpsc_executable_name = vpsc_name;
                run.magic_vpsc_box_path = vpsc_path;
                run.write_inputs_to_magic_vpsc_box;
                run.call_VPSC_executable;
                run.cleanup;
                old_run = run;

                % Set it back to normal iSave values.
                old_run.parameters.iSave = 11;
                old_run.parameters.nProcess = 1;
                old_run.parameters.processType = 0;
                old_run.parameters.processDetail(2,:) = [];

                % Get ODF from VPSC run
                outputODF(it,1) = calcDensity(old_run.texture_out{1}.orientations,'weights',old_run.texture_out{1}.weights, 'kernel', psi);
                outputODF(it,2) = calcDensity(old_run.texture_out{2}.orientations,'weights',old_run.texture_out{2}.weights, 'kernel', psi);
                cumulativeRot(it) = cumulativeRot(it-1)*rotation.byMatrix(rot_mat);

                activities.prism(it) = 0;
                activities.basal(it) = 0;
                activities.pyra(it)  = 0;
                activities.pyrca(it) = 0;
        end

    end
    %%

    % =========================================================================
    % --- SAVE and Cleanup ---
    % =========================================================================

    output_directory = 'outputs\1600_2p';
    mkdir(fullfile(MainDir,output_directory));

    activities.sum = ...
        activities.prism + ...
        activities.basal + ...
        activities.pyra  + ...
        activities.pyrca;

    %This saves the data to a MATLAB .mat file which gets loaded by the analysis script.
    % pointTrackFile(1:end-4) strips the ".csv" suffix, but if the file is .xslx
    % instead, it would need to be (1:end-5). Using a file that's been converted in
    % excel also causes other issue. Using the CSV directly from DEFORM is better.
    save([output_directory filesep pointTrackFile(1:end-4) '_Point_' num2str(pt) '.mat'], 'outputODF', 'phaseFraction', 'cumulativeRot', 'step', 'activities')

end

disp("DONE!!!!!")

toc

Analysis3;