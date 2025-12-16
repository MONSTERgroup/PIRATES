function [vel_subset,strain_inc_subset,temp_subset,time_subset,step_subset, sr_subset] = ...
    subsetDefPTR(j,time,temp,strain_inc,vel_grad,step, strain_rate)

%Select quantities of interest for current point only
        vel_subset = squeeze(vel_grad(:,j,:));
    strain_inc_subset = squeeze(strain_inc(:,j,:));
    temp_subset = squeeze(temp(:,j));
    time_subset = squeeze(time(:,j));
    step_subset = squeeze(step(:,j));
    sr_subset = squeeze(strain_rate(:,j,:));




















