warning off MATLAB:divideByZero;
warning off MATLAB:singularMatrix;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;clc
x = ([0 45 90 135 180 225 270 315]*pi/180)'; % for fitting
for i = 1:126
    xdata(:, i) = x;    
end
t = 0:0.02:0.02*(126-1);  %% 190 bins -- ingore the last bins with zeros
for i = 1:8
    tdata(i,:) = t;
end

xtdata = [xdata; tdata];
% vect=[0.258 2.430 1.155 192.2/180*pi 1.294 0.257 0.615 0.237 169.8/180*pi 0.744 224.3/180*pi];
vect=[0.258 2.430 1.155 192.2/180*pi 1.294 0.257 0.615 0.301 169.8/180*pi 0.001 224.3/180*pi];

global model_use
model_use=3;
spacetime_data = funccosnlin(vect,xtdata);

numran=100;

for i=1:numran
    i
    Amp = max(max(spacetime_data)) - min(min(spacetime_data));
    noisy_data = spacetime_data + .2*Amp*randn(size(spacetime_data));
    
    global model_use
    model_use=1;
    [spacefit_Vel,vect_Vel,r_squared_Vel,CI_Vel] = MSF_Vel_fit(noisy_data);    
    error_surf = noisy_data - spacefit_Vel;
    err_Vel = cosnlin_err(vect_Vel);
    
    model_use=2;
    [spacefit_VelAcc,vect_VelAcc,r_squared_VelAcc,CI_VelAcc] = MSF_VelAcc_fit(noisy_data, vect_Vel);
    error_surf = noisy_data - spacefit_VelAcc;    
    err_VelAcc = cosnlin_err(vect_VelAcc);
    
    model_use=3;
    [spacefit_VelAccPos,vect_VelAccPos,r_squared_VelAccPos,CI_VelAccPos] = MSF_VelAccPos_fit(noisy_data, vect_VelAcc);
    error_surf = noisy_data - spacefit_VelAccPos;
    err_VelAccPos = cosnlin_err(vect_VelAccPos);    
    
    %F-test for comparing these two model
    Ftest_1vs2(i,1)=[(err_Vel-err_VelAcc)/(length(vect_VelAcc)-length(vect_Vel))]/[err_VelAcc/(size(noisy_data,1)*size(noisy_data,2)-length(vect_Vel))];
    p_1vs2(i,1)=1-fcdf(Ftest_1vs2(i,1),length(vect_VelAcc)-length(vect_Vel),size(noisy_data,1)*size(noisy_data,2)-length(vect_Vel));    
   
    Ftest_2vs3(i,1)=[(err_VelAcc-err_VelAccPos)/(length(vect_VelAccPos)-length(vect_VelAcc))]/[err_VelAccPos/(size(spacetime_data,1)*size(spacetime_data,2)-length(vect_VelAcc))];
    p_2vs3(i,1)=1-fcdf(Ftest_2vs3(i,1),length(vect_VelAccPos)-length(vect_VelAcc),size(spacetime_data,1)*size(spacetime_data,2)-length(vect_VelAcc));    
    save Model_1vs2vs3 Ftest_1vs2 p_1vs2 Ftest_2vs3 p_2vs3
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




    
