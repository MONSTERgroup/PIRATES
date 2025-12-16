% =========================================================================
% --- ANALYSIS SCRIPT ---
% Generates figures.
% =========================================================================

% These parameters should be the same as in Main3 wherever possible.

MainDir = pwd;
if ~exist(fullfile(pwd, 'Initial Inputs'), 'dir')
    error('Cannot find initial inputs folder; confirm that the script is being called from the correct working directory.');
end

pointTrackFile = 'DEF_PTR001_1600_2p.CSV';
folder = 'outputs\1600_2p'; % This is the output_directory in Main3.m

% This is where the figure output goes. It can be the same directory or somewhere else
% entirely.
out_folder = fullfile(MainDir,folder);
mkdir(out_folder);

% List of points you want to analyze. Does not need to be the full list of points.
wantedPt = [2];

plot_video = false;

% Names of points for labeling figures. Needs to be as long as the number of points in
% the point tracking file.
pos_name = {'ND_Offset', 'Center', 'TD_Offset', 'RD_Offset'};

% Thresholds for the pole figure color bars.
MRD_Thresh_a = 4;
MRD_Thresh_b = 4;

% These parameters should be identical to their counterparts in Main3.
CSa = crystalSymmetry('622', [2.95 2.95 4.68], 'X||a', 'Y||b*', 'Z||c', 'color', 'light blue');
CSb = crystalSymmetry('432', [3.24 3.24 3.24], 'mineral', 'Titanium - Beta', 'color', 'light green');
SS = specimenSymmetry('-1');
psi = deLaValleePoussinKernel('halfwidth', 8*degree);

% =========================================================================
% --- MTEX Plotting Conventions ---
% =========================================================================

plotx2north;
plotzOutOfPlane;

% Set PF annotations to match standard viewing directions
setMTEXpref('FontSize', 40)
pfAnnotations = @(varargin) text([vector3d.X,-vector3d.Y], {'RD', 'TD'},...
    'BackgroundColor', 'w', 'FontSize', 32, 'tag', 'axesLabels', varargin{:});
setMTEXpref('pfAnnotations', pfAnnotations);

% This may need to be changed to match the rolling direction in DEFORM.

%%

% =========================================================================
% --- PLOT PHASE FRACTION vs STEP ---
% =========================================================================

for pt = wantedPt

    load([folder filesep pointTrackFile(1:end-4) '_Point_' num2str(pt) '.mat']);

    % Plot and save the phase fractions
    figure;
    set(gcf,'Color', 'w');
    plot(phaseFraction(2:end,:))
    ylim([0 1])
    xlabel('step')
    ylabel('phase fraction')
    legend('beta','alpha')
    PrettyPlotsSingle;

    export_fig([out_folder filesep 'phaseFrac_point_' pos_name{pt} '.png'], '-png');

    figure;
    plot(activities.prism ./ activities.sum)
    hold on
    plot(activities.basal ./ activities.sum)
    plot(activities.pyra  ./ activities.sum)
    plot(activities.pyrca ./ activities.sum)
    xlabel('step')
    ylabel('slip system activity (relative)')
    legend({'prism','basal','pyr <a>','pyr <c+a>'}, 'Location','best');
    PrettyPlotsSingle

    export_fig([out_folder filesep 'slipActivity_point_' pos_name{pt} '.png'], '-png');

    clear('outputODF', 'phaseFraction', 'cumulativeRot');
end
%%

% =========================================================================
% --- PLOT POLE FIGURES ---
% =========================================================================

for pt = wantedPt

    load([folder filesep pointTrackFile(1:end-4) '_Point_' num2str(pt) '.mat']);

    % Plot and save the initial textures - lab reference frame
    figure;
    plotPDF(outputODF(1,2), Miller({0,0,0,1},{1,1,-2,0},CSa),'antipodal','projection','eangle');
    plotPDF(outputODF(1,2), Miller({0,0,0,1},{1,1,-2,0},CSa),'antipodal', 'contour', 0.5:0.5:4, 'linewidth', 1.5, 'linecolor', 'black', 'ShowText', false, 'labelspacing', 300,'projection', 'eangle', 'add2all');
    set(gcf, 'units', 'pixels', 'position', [0 0 1920 1080])
    %sgtitle(['Step ' num2str(1)]);
    mtexColorMap parula
    setColorRange([0 MRD_Thresh_a])
    mtexColorbar
    export_fig([out_folder filesep 'AlphaODF_initial_point_' pos_name{pt}], '-png', '-rgb');

    figure;
    plotPDF(outputODF(1,1), Miller({1,0,0},{1,1,0},CSb),'antipodal','projection','eangle');
    plotPDF(outputODF(1,1), Miller({1,0,0},{1,1,0},CSb),'antipodal', 'contour', 0.5:0.5:4, 'linewidth', 1.5, 'linecolor', 'black', 'ShowText', false, 'labelspacing', 300,'projection', 'eangle', 'add2all');
    set(gcf, 'units', 'pixels', 'position', [0 0 1920 1080])
    %sgtitle(['Step ' num2str(1)]);
    mtexColorMap parula
    setColorRange([0 MRD_Thresh_b])
    mtexColorbar
    export_fig([out_folder filesep 'BetaODF_initial_point_' pos_name{pt}], '-png', '-rgb');


    % Plot and save the final textures - lab reference frame
    % This currently does the final textures as a separate figure for each pole. That
    % was better for making figures for the paper. You can combine them like the above
    % version by adding multple poles to Miller() in the call to plotPDF.

    figure;
    plotPDF(outputODF(end,2), Miller({0,0,0,1},CSa),'antipodal','projection','eangle');
    plotPDF(outputODF(end,2), Miller({0,0,0,1},CSa),'antipodal', 'contour', 0.5:0.5:4, 'linewidth', 1.5, 'linecolor', 'black', 'ShowText', false, 'labelspacing', 300, 'projection','eangle', 'add2all');
    set(gcf, 'units', 'pixels', 'position', [0 0 1920 1080])
    %sgtitle(['Step ' num2str(length(phaseFraction))]);
    mtexColorMap parula
    setColorRange([0 MRD_Thresh_a])
    %mtexColorbar('location','southoutside')
    export_fig([out_folder filesep 'AlphaODF_final_point_0001_' pos_name{pt}], '-png', '-rgb');

    figure;
    plotPDF(outputODF(end,2), Miller({1,1,-2,0},CSa),'antipodal','projection','eangle');
    plotPDF(outputODF(end,2), Miller({1,1,-2,0},CSa),'antipodal', 'contour', 0.5:0.5:4, 'linewidth', 1.5, 'linecolor', 'black', 'ShowText', false, 'labelspacing', 300, 'projection','eangle', 'add2all');
    set(gcf, 'units', 'pixels', 'position', [0 0 1920 1080])
    %sgtitle(['Step ' num2str(length(phaseFraction))]);
    mtexColorMap parula
    setColorRange([0 MRD_Thresh_a])
    %mtexColorbar
    export_fig([out_folder filesep 'AlphaODF_final_point_1120_' pos_name{pt}], '-png', '-rgb');

    figure;
    plotPDF(outputODF(end,1), Miller({1,0,0},CSb),'antipodal','projection','eangle');
    plotPDF(outputODF(end,1), Miller({1,0,0},CSb),'antipodal', 'contour', 0.5:0.5:4, 'linewidth', 1.5, 'linecolor', 'black', 'ShowText', false, 'labelspacing', 300, 'projection','eangle', 'add2all');
    set(gcf, 'units', 'pixels', 'position', [0 0 1920 1080])
    %sgtitle(['Step ' num2str(length(phaseFraction))]);
    mtexColorMap parula
    setColorRange([0 MRD_Thresh_b])
    %mtexColorbar
    export_fig([out_folder filesep 'BetaODF_final_point_100_' pos_name{pt}], '-png', '-rgb');

    figure;
    plotPDF(outputODF(end,1), Miller({1,1,0},CSb),'antipodal','projection','eangle');
    plotPDF(outputODF(end,1), Miller({1,1,0},CSb),'antipodal', 'contour', 0.5:0.5:4, 'linewidth', 1.5, 'linecolor', 'black', 'ShowText', false, 'labelspacing', 300, 'projection','eangle', 'add2all');
    set(gcf, 'units', 'pixels', 'position', [0 0 1920 1080])
    %sgtitle(['Step ' num2str(length(phaseFraction))]);
    mtexColorMap parula
    setColorRange([0 MRD_Thresh_b])
    %mtexColorbar
    export_fig([out_folder filesep 'BetaODF_final_point_110_' pos_name{pt}], '-png', '-rgb');

    close all

    clear('outputODF', 'phaseFraction', 'cumulativeRot');

end
%%

% =========================================================================
% --- PLOT TEXTURE EVOLUTION AS VIDEO ---
% This is slow and inelegant, but makes videos of the texture evolution for the alpha
% phase. You can duplicate this and change the proper values for the beta phase. I
% recommend commenting the entire section out if you don't need the videos.
% =========================================================================

if plot_video

    for pt = wantedPt

        load([folder filesep pointTrackFile(1:end-4) '_Point_' num2str(pt) '.mat']);

        if exist(fullfile([out_folder filesep 'AlphaODF_gif_point_' pos_name{pt} '.avi']), 'file')
            delete(fullfile([out_folder filesep 'AlphaODF_gif_point_' pos_name{pt} '.avi']));
        end

        v = VideoWriter([out_folder filesep 'AlphaODF_gif_point_' pos_name{pt} '.avi'],'Motion JPEG AVI');
        open(v);
        for k = 1:length(outputODF)
            f=figure;
            plotPDF(outputODF(k,2), Miller({0,0,0,1},{1,1,-2,0},CSa),'antipodal','projection','eangle');
            plotPDF(outputODF(k,2), Miller({0,0,0,1},{1,1,-2,0},CSa),'antipodal', 'contour', 0.5:0.5:4, 'linewidth', 1.5, 'linecolor', 'black', 'ShowText', false, 'labelspacing', 300, 'projection','eangle', 'add2all');
            %set(gcf, 'units', 'pixels', 'position', [0 0 1920 1080]) % Normally I'd have
            %this line turned on to make sure everything is the same size, but an MTEX
            %plotting quirk means that the videos turn out better (usually) if you let it
            %choose the size of the figure on its own.
            sgtitle(['Iter ' num2str(k)]);
            mtexColorMap parula
            setColorRange([0 MRD_Thresh_a])
            mtexColorbar('location','southoutside')
            frame = getframe(f);
            writeVideo(v,frame)
            close all;
        end
        close(v);
        clear('outputODF', 'phaseFraction', 'cumulativeRot');
    end

end
