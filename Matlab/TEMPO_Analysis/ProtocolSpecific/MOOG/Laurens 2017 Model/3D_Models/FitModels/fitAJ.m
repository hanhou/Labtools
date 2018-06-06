% fit Velocity-acceleration-jerk model for 3D tuning
% time (unit: s)
% spon: spontaneous firing rate
% PSTH_data: PSTH data with sliding windows
% spatial_data: 9*5 data according to spatial
% reps: the repetition number to perform for fiiting models
% 20170603LBY

function [modelFitRespon_AJ,modelFit_AJ, modelFitPara_AJ, BIC_AJ, RSquared_AJ, rss_AJ, time] = fitAJ(spon,PSTH_data,spatial_data, nBins,reps,stimOnBin,stimOffBin,aMax,aMin)

sprintf('Fitting AJ model...')

%-- initialize global using parameters

% spatial parameters
u_azi = [0 45 90 135 180 225 270 315]';
u_ele = [-90 -45 0 45 90]';
s_data = [u_ele;u_azi]; % transform to this form for fitting
% time parameters
time = (1: (stimOffBin - stimOnBin +1))' /(stimOffBin - stimOnBin +1)*1500/1000;
st_data = [u_ele;u_azi;time]; % transform to this form for fitting
% fitting initial parameters
sig = sqrt(sqrt(2))/6;
baseline = spon;
acc_max = aMax/1000; % transfer unit from ms to s
acc_min = aMin/1000; % transfer unit from ms to s
mu = (acc_max+acc_min)/2;

% PSTH data
temporal_data = squeeze(mean(mean(PSTH_data(:,:,:),1),2)); % PSTH data according to time bin
spatial_data = permute(spatial_data,[2 1]);
y_data = permute(PSTH_data, [2 1 3]); % transform to azi*ele*timebin

d_gauss_time = acc_func([mu sig],time);
d2_gauss_time = jerk_func([mu sig],time);

u1 = d_gauss_time;
p = (u1'*d2_gauss_time)/(u1'*u1);
u2 = d2_gauss_time - p*u1;

t_psth = y_data - baseline;
a_s_profile = zeros([length(u_azi), length(u_ele)]);
j_s_profile = zeros([length(u_azi), length(u_ele)]);

for j=1:length(u_ele),
    for i=1:length(u_azi),
        t_profile = squeeze(t_psth(i,j,:));
        coeff = (pinv([u1 u2])*squeeze(t_profile));
        a_s_profile(i,j) = coeff(1)-coeff(2)*p;
        j_s_profile(i,j) = coeff(2);
    end
end

% normalise time and spatial profile
a_DC = (min(a_s_profile(:))+max(a_s_profile(:)))/2;
a_A = (max(a_s_profile(:))-min(a_s_profile(:)))/2;
a_space_profile = a_s_profile-a_DC;
a_space_profile = a_space_profile/a_A;

j_DC = (min(j_s_profile(:))+max(j_s_profile(:)))/2;
j_A = (max(j_s_profile(:))-min(j_s_profile(:)))/2;
j_space_profile = j_s_profile-j_DC;
j_space_profile = j_space_profile/j_A;

%optimisation parameters for profile fits
options = optimset('Display', 'off', 'MaxIter', 5000);


%% fit acceleration spatial profile
%-- 1st, fit spatial profile1
LB = [0.001 0 -pi/2];
UB = [10 2*pi pi/2];

[~, max_idx] = max(a_space_profile(:));
[max_idx_a, max_idx_e] = ind2sub(size(a_space_profile), max_idx);

param = [0.01 u_azi(max_idx_a) u_ele(max_idx_e)];
recon_v = lsqcurvefit('cos_tuning', param,  s_data, ...
    a_space_profile(:), LB, UB, options);
a_n = recon_v(1);
a_a_0 = recon_v(2);
a_e_0 = recon_v(3);

%Initialise spatial tuning function error
ele_azi = cos_tuning(recon_v, s_data);
ele_azi = reshape(ele_azi, length(u_azi), length(u_ele));
a_space_profile2 = a_space_profile - ele_azi;

% normalise time and spatial profile
a_DC2 = (max(a_space_profile2(:)) + min(a_space_profile2(:)))/2;
a_s_A2 = (max(a_space_profile2(:)) - min(a_space_profile2(:)))/2;
a_space_profile2 = (a_space_profile2 - a_DC2)/a_s_A2;


%-- 2nd, fit spatial profile2
[~, max_idx] = max(a_space_profile2(:));
[max_idx_a, max_idx_e] = ind2sub(size(a_space_profile2), max_idx);
param = [0.01 u_azi(max_idx_a) u_ele(max_idx_e)];
recon_v = lsqcurvefit('cos_tuning', param,  s_data, ...
    a_space_profile2(:), LB, UB, options);
a_n2 = recon_v(1);

[x0, y0, z0] = sph2cart(a_a_0, a_e_0, 1);
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
a_a2_0 = ang;
a_e2_0 = recon_v(3)-a_e_0;


%-- 3rd, fit total spatial profile

% initialize fitting
% Inital fits
param = [a_n, ...   %1
    a_a_0, ... %2
    a_e_0, ... %3
    a_n2, ...  %4
    a_a2_0, ...%5
    a_e2_0, ...%6
    a_s_A2, ...%7
    a_DC2];    %8

LB = [0.001 0 -pi/2 0.001 0 elim 0 -2];
UB = [10 2*pi pi/2  10 pi 2*pi-elim 1 2];

recon_v = lsqcurvefit('d_cos_tuning', param, s_data, ...
    a_space_profile(:), LB, UB, options);

a_n    = recon_v(1);
a_a_0  = recon_v(2);
a_e_0  = recon_v(3);
a_n2   = recon_v(4);
a_a2_0 = recon_v(5);
a_e2_0 = recon_v(6);
a_s_A2 = recon_v(7);
a_DC2  = recon_v(8);

%% fit jerk spatial profile
LB = [0.001 0 -pi/2];
UB = [10 2*pi pi/2];

[~, max_idx] = max(j_space_profile(:));
[max_idx_a, max_idx_e] = ind2sub(size(j_space_profile), max_idx);

param = [0.01 u_azi(max_idx_a) u_ele(max_idx_e)];
recon_v = lsqcurvefit('cos_tuning', param,  s_data, ...
    j_space_profile(:), LB, UB, options);
j_n = recon_v(1);
j_a_0 = recon_v(2);
j_e_0 = recon_v(3);

%Initialise spatial tuning function error
ele_azi = cos_tuning(recon_v, s_data);
ele_azi = reshape(ele_azi, length(u_azi), length(u_ele));
j_space_profile2 = j_space_profile - ele_azi;

% normalise time and spatial profile
j_DC2 = (max(j_space_profile2(:)) + min(j_space_profile2(:)))/2;
j_s_A2 = (max(j_space_profile2(:)) - min(j_space_profile2(:)))/2;
j_space_profile2 = (j_space_profile2 - j_DC2)/j_s_A2;


% 2nd, fit spatial profile2
[~, max_idx] = max(j_space_profile2(:));
[max_idx_a, max_idx_e] = ind2sub(size(j_space_profile2), max_idx);
param = [0.01 u_azi(max_idx_a) u_ele(max_idx_e)];
recon_v = lsqcurvefit('cos_tuning', param,  s_data, ...
    j_space_profile2(:), LB, UB, options);
j_n2 = recon_v(1);

[x0, y0, z0] = sph2cart(j_a_0, j_e_0, 1);
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
j_a2_0 = ang;
j_e2_0 = recon_v(3)-j_e_0;


% 3rd, fit total spatial profile

% initialize fitting
% Inital fits
param = [j_n, ...   %1
    j_a_0, ... %2
    j_e_0, ... %3
    j_n2, ...  %4
    j_a2_0, ...%5
    j_e2_0, ...%6
    j_s_A2, ...%7
    j_DC2];    %8

LB = [0.001 0 -pi/2 0.001 0 elim 0 -2];
UB = [10 2*pi pi/2  10 pi 2*pi-elim 1 2];

recon_v = lsqcurvefit('d_cos_tuning', param, s_data, ...
    j_space_profile(:), LB, UB, options);

j_n    = recon_v(1);
j_a_0  = recon_v(2);
j_e_0  = recon_v(3);
j_n2   = recon_v(4);
j_a2_0 = recon_v(5);
j_e2_0 = recon_v(6);
j_s_A2 = recon_v(7);
j_DC2  = recon_v(8);
%%

a_DC = a_DC/a_A;
j_DC = j_DC/j_A;

%Fit linear parameters
A = a_A+j_A;
R_0 = baseline;
w_a = a_A/A;

%% fit AJ model
mu_t = mu;
sig_t = sig;

%Inital fits
param = [A, ...       %1
    R_0, ...     %2
    mu_t, ...    %3
    sig_t, ...   %4
    a_n, ...       %5
    a_a_0, ...     %6
    a_e_0, ...     %7
    a_n2, ...      %8
    a_a2_0, ...    %9
    a_e2_0, ...    %10
    a_s_A2, ...    %11
    a_DC+a_DC2,... %12
    j_n, ...           %13
    j_a_0, ...         %14
    j_e_0, ...         %15
    j_n2, ...          %16
    j_a2_0, ...        %17
    j_e2_0, ...        %18
    j_s_A2, ...        %19
    j_DC2+j_DC, ...%20
    w_a];               %21

init_param = zeros(reps+1, length(param));
init_param(1,:) = param;

LB = [0.25*A, ...`  %1  A
    0, ...          %2  R_0
    mu-0.1, ...       %3  mu_t
    0.5*sig, ...    %4  sig_t
    0.001, ...      %5  n
    0, ...          %6  a_0
    -pi/2, ...      %7  e_0
    0.001, ...      %8  n2
    0, ...          %9  a2_0
    elim, ...       %10 e2_0
    0, ...          %11 s_A2
    -2,...          %12 v_DC
    0.001, ...      %13 a_n
    0, ...          %14 j_a_0
    -pi/2, ...      %15 j_e_0
    0.001, ...      %16 j_n2
    0, ...          %17 j_a2_0
    elim, ...       %18 j_e2_0
    0, ...          %19 j_s_A2
    -2, ...         %20 j_DC
    0];             %21   w_a

UB = [4*A, ...      %1  A
    300, ...        %2  R_0
    mu+0.1, ...      %3  mu_t
    2*sig, ...      %4  sig_t
    10, ...         %5  n
    2*pi, ...       %6  a_0
    pi/2, ...       %7  e_0
    10, ...         %8  n2
    pi, ...         %9  a2_0
    2*pi-elim, ...  %10 e2_0
    1, ...          %11 s_A2
    2,...          %12 v_DC
    10, ...        %13 j_n
    2*pi, ...      %14 j_a_0
    pi/2, ...      %15 j_e_0
    10, ...        %16 j_n2
    pi, ...        %17 j_a2_0
    2*pi-elim, ... %18 j_e2_0
    1, ...         %19 j_s_A2
    2, ...         %20 j_DC
    1];            %21    w_a


rand_rss = zeros(reps+1,1);
rand_param = zeros(reps+1, length(param));
rand_jac = zeros(reps+1, length(param), length(param));

[rand_param(1,:),rand_rss(1),~,~,~,~,temp_jac] = lsqcurvefit('AJ_Model', ...
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
    
    [rand_param(ii,:),rand_rss(ii),~,~,~,~,temp_jac] = lsqcurvefit('AJ_Model', ...
        seed_param, st_data, y_data, LB, UB, options);
    rand_jac(ii,:,:) = full(temp_jac)'*full(temp_jac);
    
    if rand_rss(ii) < min_rss
        min_rss = rand_rss(ii);
        min_param = rand_param(ii,:);
    end
    
    
end

% find the best fit parameters according to rss
[~,min_inx] = min(rand_rss);
modelFitPara_AJ = rand_param(min_inx,:);
rss_AJ = rand_rss(min_inx);
jac_AJ = rand_jac(min_inx,:,:);

% calculate the final model fitting values
respon = AJ_Model(modelFitPara_AJ,st_data);
modelFitRespon_AJ = respon;

modelFit_AJ.A = AO_Model(modelFitPara_AJ(1:12),st_data);
modelFit_AJ.J = JO_Model(modelFitPara_AJ([1:4 13:20]),st_data);
%% analysis
data_num = 26*nBins;
para_num = 10;
BIC_AJ = BIC_fit(data_num,rss_AJ,para_num);
TSS = sum((PSTH_data(:) - mean(PSTH_data(:))).^2);
RSquared_AJ = 1 - rss_AJ/TSS;

end
