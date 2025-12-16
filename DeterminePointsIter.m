% Function to determine the number of points and number of steps per point
% in a deformation file. 

function [pt,iter] = DeterminePointsIter(fname)

data = xlsread(fname);

pt = unique(data(:,1),'stable');
iter = round(length(data(:,2))/length(pt));

end