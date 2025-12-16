% Function to read in relevant variables from DEFORM *.csv file
% V.M. Miller July 2016
% DEFORM version 11.1 PTR file compatible. 

function [vel_grad,strain_inc, temp, time, pos,step, strain_rate] = readDEFPTR3d(fname,pt,iter)

data = xlsread(fname);

strain_incx = reshape(data(:,106),iter,pt);
strain_incy = reshape(data(:,107),iter,pt);
strain_incz = reshape(data(:,108),iter,pt);

temp = reshape(data(:,15),iter,pt);

time = reshape(data(:,4),iter,pt);

step = reshape(data(:,2),iter,pt);

% velocity gradient is 123 to 131 <- this is right 
% deformation gradient is 132 to 140 <- this is wrong
vg1x = reshape(data(:,123),iter,pt);
vg1y = reshape(data(:,124),iter,pt);
vg1z = reshape(data(:,125),iter,pt);
vg2x = reshape(data(:,126),iter,pt);
vg2y = reshape(data(:,127),iter,pt);
vg2z = reshape(data(:,128),iter,pt);
vg3x = reshape(data(:,129),iter,pt);
vg3y = reshape(data(:,130),iter,pt);
vg3z = reshape(data(:,131),iter,pt);

srxx = reshape(data(:,37),iter,pt);
srxy = reshape(data(:,40),iter,pt);
srxz = reshape(data(:,42),iter,pt);
sryx = reshape(data(:,40),iter,pt);
sryy = reshape(data(:,38),iter,pt);
sryz = reshape(data(:,41),iter,pt);
srzx = reshape(data(:,42),iter,pt);
srzy = reshape(data(:,41),iter,pt);
srzz = reshape(data(:,39),iter,pt);


x = reshape(data(:,8),iter,pt);
y = reshape(data(:,9),iter,pt);
z = reshape(data(:,10),iter,pt);

pos = zeros(iter,pt,3);
pos(:,:,1) = x;
pos(:,:,2) = y;
pos(:,:,3) = z;

strain_inc = zeros(iter,pt,3);
strain_inc(:,:,1) = strain_incx;
strain_inc(:,:,2) = strain_incy;
strain_inc(:,:,3) = strain_incz;

vel_grad = zeros(iter,pt,9);
vel_grad(:,:,1)= vg1x; %xx
vel_grad(:,:,2)= vg1y; %xy
vel_grad(:,:,3)= vg1z; %xz
vel_grad(:,:,4)= vg2x; %yx
vel_grad(:,:,5)= vg2y; %yy
vel_grad(:,:,6)= vg2z; %yz
vel_grad(:,:,7)= vg3x; %zx
vel_grad(:,:,8)= vg3y; %zy
vel_grad(:,:,9)= vg3z; %zz

strain_rate = zeros(iter,pt,9);
strain_rate(:,:,1)= srxx; %xx
strain_rate(:,:,2)= srxy; %xy
strain_rate(:,:,3)= srxz; %xz
strain_rate(:,:,4)= sryx; %yx
strain_rate(:,:,5)= sryy; %yy
strain_rate(:,:,6)= sryz; %yz
strain_rate(:,:,7)= srzx; %zx
strain_rate(:,:,8)= srzy; %zy
strain_rate(:,:,9)= srzz; %zz

end
