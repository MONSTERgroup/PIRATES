function [curr_vel,curr_inc,delta_T, curr_sr] = getIncrVals(it,vel_subset,strain_inc_subset,temp_subset, sr_subset,varargin)

dimension = 3;

        % Write deformation file
        if dimension == 2
            curr_rate = rate_subset(it,:); %
        elseif dimension == 3
            curr_vel = vel_subset(it,:);
            curr_sr = sr_subset(it,:);
        end

        curr_inc = (strain_inc_subset(it,:)-strain_inc_subset(it-1,:));
        delta_T = temp_subset(it)-temp_subset(it-1);