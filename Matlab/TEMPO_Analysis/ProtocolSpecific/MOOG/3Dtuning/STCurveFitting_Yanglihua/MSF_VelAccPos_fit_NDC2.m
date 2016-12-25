function [space_fit,vect_calculated,r_squared,Current_ci,cor] = MSF_VelAccPos_fit_NDC2(spacetime_data, vect_VelAcc,timestep,allow_negative)
% MSFfit fits a Modified Sinusoid Function (MSF) to data using 'fmincon', with parameter bounds
%	It calls MSFerr.m to comput the error for each set of params.  The function 
%   evaluated is given by MSFfunc.m

clear global rawdata xdata tdata
global rawdata xdata tdata 
rawdata=spacetime_data;

x = ([0 45 90 135 180 225 270 315 360]*pi/180)'; % for fitting
for i = 1:size(spacetime_data,2)
    xdata(:, i) = x;    
end

t=0:timestep:timestep*(size(spacetime_data,2)-1);
for i = 1:size(spacetime_data,1)
    tdata(i,:) = t;
end
xtdata = [xdata; tdata];  

max_val=max(max(rawdata));
min_val=min(min(rawdata));
max_indx=find(rawdata==max_val);
min_indx=find(rawdata==min_val);
maxazi = find( mean(rawdata,2) == max(mean(rawdata,2)) );
maxazi_data = rawdata(maxazi, :);

global model_use
model_use=3;

min_err = 9999999999999999.99;
if allow_negative
    LB = [0 -1.5*(max_val-min_val) 0.0001 (0*45*pi/180) 0 0 0 (0*45*pi/180) (0*45*pi/180) 0];   % lower bounds
    UB = [1.35*max_val 1.5*(max_val-min_val) 10 (8*45*pi/180) 2.5 6 1 (8*45*pi/180) (8*45*pi/180) 1];   % upper bounds
else
    LB = [0 0 0.0001 (0*45*pi/180) 0 0 0 (0*45*pi/180) (0*45*pi/180) 0];   % lower bounds
    UB = [1.35*max_val 1.5*(max_val-min_val) 10 (8*45*pi/180) 2.5 6 1 (8*45*pi/180) (8*45*pi/180) 1];   % upper bounds
end

N = 10;
vect = [vect_VelAcc 0 0];
vect_temp7 = LB(7) : (UB(7)-LB(7))/(N-1) : UB(7);
vect_temp9 = LB(9) : (UB(9)-LB(9))/(N-1) : UB(9);
vect_temp10 = LB(10) : (UB(10)-LB(10))/(N-1) : UB(10);

                  
for jj=1:N
    vect_t9 = vect_temp9(jj);%     
    for nn=1:N  
        vect_t10 = vect_temp10(nn);
        vect_temp = [vect(1) vect(2) vect(3) vect(4) vect(5) vect(6) vect(7) vect(8) vect_t9 vect_t10];
        error = cosnlin_err_NDC2(vect_temp);
        if (error < min_err) 
            nn
            vect_t9min = vect_t9;
            vect_t10min = vect_t10;                   
            min_err = error;
        end
    end
end

if error<=min_err  
    vect(9) = vect_t9min;
    vect(10) = vect_t10min;    
else
    vect = [vect_VelAcc 0 0];
end

min_err = cosnlin_err_NDC2(vect);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ydata = rawdata;          
wiggle = 0.9;                  
N_reps = 40;
options = optimset('Display', 'off', 'MaxIter', 5000, 'LevenbergMarquardt', 'on'); % 'LevenbergMarquardt', 'on'); % 'Tolx',1e-4,
A = []; b = []; Aeq = []; beq = []; nonlcon = [];

for j=1:N_reps 
%     j
    rand_factor = rand(size(vect)) * wiggle + (1-wiggle/2); % ranges from 1-wiggle/2 -> 1 + wiggle/2
    temp_vect = vect .* rand_factor;
    [testpars{j},resnorm{j},residual{j},exitflag{j},output{j},lambda{j},jacobian{j}] = lsqcurvefit('funccosnlin_NDC2', temp_vect, xtdata, ydata, LB, UB, options);
    err_pars(j) = cosnlin_err_NDC2(testpars{j});    
    clear rand_factor temp_vect;
end

[min_err min_indx] = min(err_pars);
clear vect_calculated;vect_calculated=testpars{min_indx};
space_fit = (funccosnlin_NDC2(vect_calculated, xtdata));
Current_err = cosnlin_err_NDC2(vect_calculated);
err_total = sum( sum(( ydata - mean(mean(ydata)) ) .^2) );
r_squared = (1 - ((Current_err)^2 / err_total));
Current_ci = nlparci(vect_calculated,residual{min_indx},jacobian{min_indx});
%%%%%%%%%%%%%%
%calculate variance cov matrix and correlation matrix of parameters
jac=full(jacobian{min_indx});
xtx=jac'*jac;
xtxinv=inv(xtx);
varinf=diag(xtxinv);
cor=xtxinv./sqrt(varinf*varinf');

return;
