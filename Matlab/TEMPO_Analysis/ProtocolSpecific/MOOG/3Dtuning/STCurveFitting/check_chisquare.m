clear;clc
load ChiSquareData.mat

x = ([0 45 90 135 180 225 270 315 360]*pi/180)'; % for fitting
for i = 1:13
    xdata(:, i) = x;    
end
t = 0:0.25:2.5;  %% 190 bins -- ingore the last bins with zeros
for i = 1:9
    tdata(i,:) = t;
end
xtdata = [xdata; tdata];

global model_use
model_use=1;
vect_Vel=[12.1 25.9 0.9 197.2*pi/180 1.039 0.208 0.533];
spacetime_data = funccosnlin(vect_Vel,xtdata);

% Do chi-square goodness of fit test
[spacefit_Vel,newvect] = funccosnlin(vect_Vel,xtdata); 
spacefit_Vel(spacefit_Vel<0)=0;
clear DataDiff;DataDiff=spacetime_data_n-sqrt(spacefit_Vel);