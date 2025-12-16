% Function to segment the point tracking file
% I think we need a tolerance to the coords equal condition to fix minor
% issue. Requires testing.

function [pt,seg_mat] = SegmentPointTrack(fname, pos)

data = xlsread(fname);
pt = unique(data(:,1),'stable');

% find the number of rows that corresponds to each point
nstep = round(length(data(:,2))/length(pt));


% Let's build a segment matrix
seg_mat = cell(nstep,1);


%Segmet File Values:
% 0 means nothing happens
% 1 means hold (update the phase fractions)
% 2 means deform (update phase fractions and run VPSC)
% 3 means rotation (update texture with rotation)

i=1;

while i <= nstep % i is row number in the data file, limited to the rows for the first point
    
    
    
    if i == 1  % the initial condition
        seg_mat{i,1} = 'nothing';
        i=i+1;
        continue
    end
    
    if data(i,8)==data(i-1,8) && ...
            data(i,9)==data(i-1,9) && ...
            data(i,10)==data(i-1,10) ...
            && data(i,2)<0 % position didn't change and step # is neg
        %just a re-mesh or boundary step
        seg_mat{i,1} = 'nothing';
        i=i+1;
        continue
        
    elseif (data(i,8)~=data(i-1,8) || ...
            data(i,9)~=data(i-1,9) || ...
            data(i,10)~=data(i-1,10)) ...
            && data(i,4)>data(i-1,4) %if the position changes and the time is running, it must be a deformation step.
        seg_mat{i,1} = 'deformation';
        i = i+1;
        continue
        
    elseif data(i,4)>data(i-1,4) %Time runs but stroke is zero
        % hold
        seg_mat{i,1} = 'hold';
        i = i+1;
        continue
        
    elseif data(i,4)==data(i-1,4) && ...
            data(i,2)<0 && ...
            (data(i,8)~=data(i-1,8) || ...
            data(i,9)~=data(i-1,9) || ...
            data(i,10)~=data(i-1,10))
        %step is negative, time doesn't run, and coords change
        %rotation
        
        seg_mat{i,1} = 'rotation';
        i = i+1;
        continue
        
    else
        seg_mat{i,1} = 'nothing';
        i = i+1;
    end
    
    % Vision: take this matrix, delete the zeros (from both this and the data matrix,
    % then count the number of  segments. Outside this function, create
    % a loop over segments that  grabs rows based on this matrix.
    % subroutines as necessary.
    
end