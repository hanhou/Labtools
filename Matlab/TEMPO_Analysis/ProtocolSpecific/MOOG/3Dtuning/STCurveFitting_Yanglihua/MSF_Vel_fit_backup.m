function [space_fit,vect_calculated r_squared,Current_ci,cor] = MSF_Vel_fit(spacetime_data,timestep,allow_negative)
% MSFfit fits a Modified Sinusoid Function (MSF) to data using 'fmincon', with parameter bounds
%	It calls MSFerr.m to comput the error for each set of params.  The function 
%   evaluated is given by MSFfunc.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear global rawdata xdata tdata
global rawdata xdata tdata %model_use
rawdata=spacetime_data;

% x = ([0 45 90 135 180 225 270 315]*pi/180)'; % for fitting
x = ([0 45 90 135 180 225 270 315 360]*pi/180)'; % for fitting
for i = 1:size(spacetime_data,2)%101
    xdata(:, i) = x;    
end
% t=0:0.02:0.02*(size(spacetime_data,2)-1);%t = 0:0.02:2;  %% 190 bins -- ingore the last bins with zeros
t=0:timestep:timestep*(size(spacetime_data,2)-1);% t=0:0.05:2.5;
for i = 1:size(spacetime_data,1)
    tdata(i,:) = t;
end
xtdata = [xdata; tdata]; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first, generate some initial parameter guesses
max_val=max(max(rawdata));
min_val=min(min(rawdata));
% max_indx=find(rawdata==max_val);max_indx=max_indx(1,1);
% min_indx=find(rawdata==min_val);min_indx=min_indx(1,1);

clear space_mean; space_mean=nanmean(rawdata');
N=sum(space_mean);
mean=sum(space_mean.*x')/N;
var=sum(((x'-mean).^2.*space_mean))/(N-1);
sigma_azi=sqrt(var);
azi_max_indx=find(space_mean==max(space_mean)); azi_max_indx=azi_max_indx(1,1);
azi_min_indx=find(space_mean==min(space_mean)); azi_min_indx=azi_min_indx(1,1);
mu_azi=x(azi_max_indx);

clear t_mean; t_mean=nanmean(rawdata);
N=sum(t_mean);
mean=sum(t_mean.*t)/N;
var=sum(((t-mean).^2.*t_mean))/(N-1);
sigma_t=sqrt(var);
t_max_indx=find(t_mean==max(t_mean)); t_max_indx=t_max_indx(1,1);
t_min_indx=find(t_mean==min(t_mean)); t_min_indx=t_min_indx(1,1);
mu_t=t(t_max_indx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vect=[];%[Ro  Amplitude n mu1  mu2  sigma_t DC2]
% vect(1)= nanmean(nanmean(rawdata));%
vect(1)=min_val; % DC overall
vect(2)=(max_val-min_val);  %amplitude (peak-trough modulation);
vect(3)=0.01; %nonlinearity
vect(4)=mu_azi;%azimuth of max response
vect(5)=mu_t; % time of max response
vect(6)=sigma_t;%vect(6)=6; %Sigma_t
vect(7)= nanmean(nanmean(rawdata))/(max_val-min_val);%DC2

global model_use
model_use=1;
% space_fit = funccosnlin(vect,xtdata);
% min_err = cosnlin_err(vect);
% figure;contourf(space_fit');colorbar;
% for i=1:8;
%     figure(1);subplot(2,4,i);plot(t,space_fit(i,:));
% end
% 
% clear temp_min_val; temp_min_val=min(min(space_fit));
% if temp_min_val<0
%     space_fit=space_fit-min_val;
%     vect(1)=vect(1)-temp_min_val;
% else
% end
% space_fit_2 = funccosnlin(vect,xtdata);
% figure; contourf(space_fit_2');colorbar;

N=20;
%Starting here, search for better starting values of vect(3), vect(6) and vect(7)
vect3temp=vect(3);vect3range=10; 
vect6temp=vect(6);vect6range=(max(t)-min(t));
vect7temp=vect(7);vect7range=1;
min_err = 9999999999999999.99;
for i=1:N
    vect(3)=vect3range/(N-1)*(i-1);
    for j=1:N
        vect(6)=min(t)+(i-1)*vect6range/(N-1);
        for k=1:N
            vect(7)=vect7range/(N-1)*(i-1);
            error=cosnlin_err(vect);
            if (error<min_err)
                vect3_min=vect(3);
                vect6_min=vect(6);
                vect7_min=vect(7);
                min_err=error;
            end            
        end
    end
end
vect(3)=vect3_min;
vect(6)=vect6_min;
vect(7)=vect7_min;
pos_err=min_err;

% space_fit_2 = funccosnlin(vect,xtdata);
% figure; contourf(space_fit_2');colorbar;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (allow_negative)
    %now check if error would be smaller for negative -going response. If
    %so, start there
    vectn(1)=max_val; % DC overall
    vectn(2)=-(max_val-min_val);  %amplitude (peak-trough modulation);
    vectn(3)=1; %nonlinearity
    vectn(4)=x(azi_min_indx);%azimuth of max response
    vectn(5)=t(t_min_indx); % time of max response
    vectn(6)=sigma_t; %Sigma_t
    vectn(7)= nanmean(nanmean(rawdata))/(max_val-min_val);%DC2
    
    %starting here, search for better starting values of vect(3), vect(6)
    vect3temp=vectn(3);vect3range=10; 
    vect6temp=vectn(6);vect6range=(max(t)-min(t));
%     min_err = 9999999999999999.99;
    for i=1:N
        vectn(3)=vect3range/N*i;
        for j=1:N
            vectn(6)=min(t)+(i-1)*vect6range/(N-1);
            for k=1:N
                vectn(7)=vect7range/(N-1)*(i-1);
                error=cosnlin_err(vectn);
                if (error<min_err)
                    vect3_min=vectn(3);
                    vect6_min=vectn(6);
                    vect7_min=vectn(7);
                    min_err=error;
                end
            end
        end
    end
    vectn(3)=vect3_min;
    vectn(6)=vect6_min;
    neg_err=min_err;  
    if neg_err<pos_err
        vect=vectn;
    end
end

if allow_negative
    LB = [0 -1.5*(max_val-min_val) 0.001 (0*45*pi/180) 0 0 -1];   %     lower bounds
    UB = [1.35*max_val 1.5*(max_val-min_val) 10 (8*45*pi/180) 2.5 6 1];   % upper bounds
else
    LB = [0 0 0.001 (0*45*pi/180) 0 0 -1];   % lower bounds
    UB = [1.35*max_val 1.5*(max_val-min_val) 10 (8*45*pi/180) 2.5 6 1];   % upper bounds
end


ydata = rawdata;         %% your results to pass to lsqcurvefit
xtdata = [xdata; tdata];          %% for lsqcurvefit, need to arrange this way since we are comparing surface
wiggle = 0.5;   % wiggle = 0.2;                  
N_reps = 100;% N_reps = 20;

options = optimset('Display', 'off', 'MaxIter', 5000, 'LevenbergMarquardt', 'on'); % 'LevenbergMarquardt', 'on'); % 'Tolx',1e-4,
A = []; b = []; Aeq = []; beq = []; nonlcon = [];

%% lsqcurvefit with jitter added to initial guesses
% min_err = 9999999999999999.99;
for j=1:N_reps
%     j
    rand_factor = rand(size(vect)) * wiggle + (1-wiggle/2); % ranges from 1-wiggle/2 -> 1 + wiggle/2
    temp_vect = vect .* rand_factor;
    [testpars{j},resnorm{j},residual{j},exitflag{j},output{j},lambda{j},jacobian{j}] = lsqcurvefit('funccosnlin', temp_vect, xtdata, ydata, LB, UB, options);
    err_pars(j) = cosnlin_err(testpars{j}); 
    clear rand_factor, temp_vect;    
end

[min_err min_indx] = min(err_pars);
clear vect_calculated;vect_calculated=testpars{min_indx};%vect_calculated = testpars_min;
space_fit = (funccosnlin(vect_calculated, xtdata));
Current_err = cosnlin_err(vect_calculated);
err_total = sum( sum(( ydata - nanmean(nanmean(ydata)) ) .^2) );
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