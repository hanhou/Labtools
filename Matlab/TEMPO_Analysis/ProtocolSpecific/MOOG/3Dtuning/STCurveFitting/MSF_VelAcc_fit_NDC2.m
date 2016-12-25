function [space_fit,vect_calculated,r_squared,Current_ci,cor] = MSF_VelAcc_fit_NDC2(spacetime_data, vect_Vel,timestep,allow_negative)
% MSFfit fits a Modified Sinusoid Function (MSF) to data using 'fmincon', with parameter bounds
%	It calls MSFerr.m to comput the error for each set of params.  The function 
%   evaluated is given by MSFfunc.m

clear global rawdata xdata tdata
global rawdata xdata tdata 
rawdata=spacetime_data;

x = ([0 45 90 135 180 225 270 315 360]*pi/180)'; % for fitting
for i = 1:size(spacetime_data,2)%101
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

global model_use
model_use=2;

min_err = 9999999999999999.99;
if allow_negative
    LB = [0 -5*(max_val-min_val) 0.0001 (0*45*pi/180) 0 0 0 (0*45*pi/180)];   % lower bounds
    UB = [1.35*max_val 1.5*(max_val-min_val) 10 (8*45*pi/180) 2.5 6 1 (8*45*pi/180)];   % upper bounds
else
    LB = [0 0 0.0001 (0*45*pi/180) 0 0 0 (0*45*pi/180)];   % lower bounds
    UB = [1.35*max_val 1.5*(max_val-min_val) 10 (8*45*pi/180) 2.5 6 1 (8*45*pi/180)];   % upper bounds
end

N = 20;
vect = [vect_Vel 0 0];
vect_temp3 = LB(3) : (UB(3)-LB(3))/(N-1) : UB(3);
vect_temp7 = LB(7) : (UB(7)-LB(7))/(N-1) : UB(7);
vect_temp8 = LB(8) : (UB(8)-LB(8))/(N-1) : UB(8);

for i = 1:N
    vect_t3 = vect_temp3(i); 
    for j = 1:N
        vect_t7 = vect_temp7(j); 
        for k = 1:N
            vect_t8 = vect_temp8(k);            
            vect_temp = [vect(1) vect(2) vect_t3 vect(4) vect(5) vect(6) vect_t7 vect_t8];
            error = cosnlin_err_NDC2(vect_temp);
            if (error < min_err)
                vect_t3min = vect_t3;
                vect_t7min = vect_t7;
                vect_t8min = vect_t8;
                min_err = error;
            end
        end
    end
end

vect(3) = vect_t3min;
vect(7) = vect_t7min;
vect(8) = vect_t8min;
min_err = cosnlin_err_NDC2(vect);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ydata = rawdata;
wiggle = 0.9;                  %wiggle = 0.9;                  
N_reps = 100;

options = optimset('Display', 'off', 'MaxIter', 5000, 'LevenbergMarquardt', 'on'); % 'LevenbergMarquardt', 'on'); % 'Tolx',1e-4,
A = []; b = []; Aeq = []; beq = []; nonlcon = [];
err_pars=[];
for j=1:N_reps
    j;
    rand_factor = rand(size(vect)) * wiggle + (1-wiggle/2); % ranges from 1-wiggle/2 -> 1 + wiggle/2
    clear temp_vect;temp_vect = vect .* rand_factor;
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