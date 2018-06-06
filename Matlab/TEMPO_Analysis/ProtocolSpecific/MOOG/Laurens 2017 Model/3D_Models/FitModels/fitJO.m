% fit Jerk-only model for 3D tuning
% time (unit: s)
% spon: spontaneous firing rate
% PSTH_data: PSTH data with sliding windows
% spatial_data: 9*5 data according to spatial
% reps: the repetition number to perform for fiiting models
% 20170603LBY

function [modelFitRespon_JO, modelFit_JO, modelFitPara_JO, BIC_JO, RSquared_JO, rss_JO, time] = fitJO(spon,PSTH_data,spatial_data, nBins,reps,stimOnBin,stimOffBin,aMax,aMin)

sprintf('Fitting JO model...')

%-- initialize global using parameters

% spatial parameters
u_azi = [0 45 90 135 180 225 270 315]';
u_ele = [-90 -45 0 45 90]';
s_data = [u_ele;u_azi]; % transform to this form for fitting
% time parameters
time = (1: (stimOffBin - stimOnBin +1))' /(stimOffBin - stimOnBin +1)*1500/1000;
st_data = [u_ele;u_azi;time]; % transform to this form for fitting
% fitting initial parameters
mu = 0.91;
sig = sqrt(sqrt(2))/6;
baseline = spon;
acc_max = aMax/1000;
acc_min = aMin/1000;

% PSTH data
temporal_data = squeeze(mean(mean(PSTH_data(:,:,:),1),2)); % PSTH data according to time bin
spatial_data = permute(spatial_data,[2 1]);
y_data = permute(PSTH_data, [2 1 3]); % transform to azi*ele*timebin

% check if the response is excitory or inhibitory
d2_gauss_time = jerk_func([mu sig],time);
corrcoeff_vel = xcorr(d2_gauss_time, temporal_data, 'coeff');
[~,max_val] = max(abs(corrcoeff_vel));
if corrcoeff_vel(max_val) < 0,
    temporal_data = -temporal_data;
    spatial_data = -spatial_data;
end

% normalise time and spatial profile
t_A = max(temporal_data) - min(temporal_data);
s_DC = (max(spatial_data(:)) + min(spatial_data(:)))/2;
s_A = (max(spatial_data(:)) - min(spatial_data(:)))/2;
temporal_data = temporal_data/t_A;
spatial_data = (spatial_data - s_DC)/s_A;

%optimisation parameters for profile fits
options = optimset('Display', 'off', 'MaxIter', 5000);

%% initial fitting before fitting model
%-- 1st, fit time profile
mu_t = mu;
sig_t = sig;

LB = [acc_max 0.5*sig];
UB = [acc_min 2*sig];
recon_t = lsqcurvefit('vel_func', [mu_t sig_t], ...
    time, temporal_data, LB, UB, options);
obj.mu_t = recon_t(1);
obj.sig_t = recon_t(2);

%-- 2nd, fit spatial profile1

LB = [0.001 0 -pi/2];
UB = [10 2*pi pi/2];

[~, max_idx] = max(spatial_data(:));
[max_idx_a, max_idx_e] = ind2sub(size(spatial_data), max_idx);

param = [0.01 u_azi(max_idx_a) u_ele(max_idx_e)];
recon_v = lsqcurvefit('cos_tuning', param,  s_data, ...
    spatial_data(:), LB, UB, options);
n = recon_v(1);
a_0 = recon_v(2);
e_0 = recon_v(3);

%Initialise spatial tuning function error
ele_azi = cos_tuning(recon_v, s_data);
ele_azi = reshape(ele_azi, length(u_azi), length(u_ele));
spatial_rate2 = spatial_data - ele_azi;

% normalise time and spatial profile
s_DC2 = (max(spatial_rate2(:)) + min(spatial_rate2(:)))/2;
s_A2 = (max(spatial_rate2(:)) - min(spatial_rate2(:)))/2;
spatial_rate2 = (spatial_rate2 - s_DC2)/s_A2;


%-- 3rd, fit spatial profile2
[~, max_idx] = max(spatial_rate2(:));
[max_idx_a, max_idx_e] = ind2sub(size(spatial_rate2), max_idx);
param = [0.01 u_azi(max_idx_a) u_ele(max_idx_e)];
recon_v = lsqcurvefit('cos_tuning', param,  s_data, ...
    spatial_rate2(:), LB, UB, options);
n2 = recon_v(1);

[x0, y0, z0] = sph2cart(a_0, e_0, 1);
[x1, y1, z1] = sph2cart(recon_v(2), recon_v(3), 1);
v1 = [1 0 0] - [x0 y0 z0];
v2 = [x1 y1 z1] - [x0 y0 z0];
ang = acos(v1*v2'/(norm(v1)*norm(v2)));

elim = pi/3;
if ang < elim,
    ang = elim;
end
if ang > 2*pi-elim,
    ang = 2*pi-elim;
end
a2_0 = ang;
e2_0 = recon_v(3)-e_0;





%-- 4th, fit total spatial profile

% initialize fitting
% Inital fits
param = [n, ...   %1
    a_0, ... %2
    e_0, ... %3
    n2, ...  %4
    a2_0, ...%5
    e2_0, ...%6
    s_A2, ...%7
    s_DC2];    %8

LB = [0.001 0 -pi/2 0.001 0 elim 0 -2];
UB = [10 2*pi pi/2  10 pi 2*pi-elim 1 2];

recon_v = lsqcurvefit('d_cos_tuning', param, s_data, ...
    spatial_data(:), LB, UB, options);

n    = recon_v(1);
a_0  = recon_v(2);
e_0  = recon_v(3);
n2   = recon_v(4);
a2_0 = recon_v(5);
e2_0 = recon_v(6);
s_A2 = recon_v(7);
s_DC2  = recon_v(8);

%Fit linear parameters
A = t_A*s_A;
R_0 = baseline;
s_DC = s_DC/s_A;

%% fit JO model

%Inital fits
param = [A, ...       %1
    R_0, ...     %2
    mu_t, ...    %3
    sig_t, ...   %4
    n, ...       %5
    a_0, ...     %6
    e_0, ...     %7
    n2, ...      %8
    a2_0, ...    %9
    e2_0, ...    %10
    s_A2, ...    %11
    s_DC+s_DC2]; %12


init_param = zeros(reps+1, length(param));
init_param(1,:) = param;

LB = [0.25*A, ...`%1  A
    0, ...          %2  R_0
    acc_max, ...      %3  mu_t
    0.5*sig, ...%4  sig_t
    0.001, ...      %5  n
    0, ...          %6  a_0
    -pi/2, ...      %7  e_0
    0.001, ...      %8  n2
    0, ...          %9  a2_0
    elim, ...       %10 e2_0
    0, ...          %11 s_A2
    -2]; ...        %12 DC
    
UB = [4*A, ...    %1  A
    300, ...        %2  R_0
    acc_min, ...        %3  mu_t
    2*sig, ... %4  sig_t
    10, ...         %5  n
    2*pi, ...       %6  a_0
    pi/2, ...       %7  e_0
    10, ...         %8  n2
    pi, ...         %9  a2_0
    2*pi-elim, ...  %10 e2_0
    1, ...          %11 s_A2
    2];             %12 DC

rand_rss = zeros(reps+1,1);
rand_param = zeros(reps+1, length(param));
rand_jac = zeros(reps+1, length(param), length(param));

[rand_param(1,:),rand_rss(1),~,~,~,~,temp_jac] = lsqcurvefit('JO_Model', ...
    init_param(1,:), st_data, y_data, LB, UB, options);
rand_jac(1,:,:) = full(temp_jac)'*full(temp_jac);
min_param = rand_param(1,:);
min_rss = rand_rss(1);
err_range =  0.1*(UB - LB);

% fitting the models
for ii = 2:(reps + 1)
    
    init_param(ii,:) = min_param;
    UB_param = min_param+err_range;
    LB_param = min_param-err_range;
    UB_param(UB < UB_param) = UB(UB < UB_param);
    LB_param(LB > LB_param) = LB(LB > LB_param);
    seed_param  = unifrnd(LB_param, UB_param);
    
    [rand_param(ii,:),rand_rss(ii),~,~,~,~,temp_jac] = lsqcurvefit('JO_Model', ...
        seed_param, st_data, y_data, LB, UB, options);
    rand_jac(ii,:,:) = full(temp_jac)'*full(temp_jac);
    
    if rand_rss(ii) < min_rss
        min_rss = rand_rss(ii);
        min_param = rand_param(ii,:);
    end
    
    
end

% find the best fit parameters according to rss
[~,min_inx] = min(rand_rss);
modelFitPara_JO = rand_param(min_inx,:);
rss_JO = rand_rss(min_inx);
jac_JO = rand_jac(min_inx,:,:);

% calculate the final model fitting values
respon = JO_Model(modelFitPara_JO,st_data);
modelFitRespon_JO = respon;

modelFit_JO = [];
%% analysis
data_num = 26*nBins;
para_num = 6;
BIC_JO = BIC_fit(data_num,rss_JO,para_num);
TSS = sum((PSTH_data(:) - mean(PSTH_data(:))).^2);
RSquared_JO = 1 - rss_JO/TSS;

end
