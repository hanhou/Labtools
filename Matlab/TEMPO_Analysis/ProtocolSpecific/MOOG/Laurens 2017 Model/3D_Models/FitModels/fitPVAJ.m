% fit Position-Velocity-acceleration-jerk model for 3D tuning
% time (unit: s)
% spon: spontaneous firing rate
% PSTH_data: PSTH data with sliding windows
% spatial_data: 9*5 data according to spatial
% reps: the repetition number to perform for fiiting models
% 20170603LBY

function [modelFitRespon_PVAJ, modelFit_PVAJ, modelFitPara_PVAJ, BIC_PVAJ, RSquared_PVAJ, rss_PVAJ, time] = fitPVAJ(spon,PSTH_data,spatial_data, nBins,reps,tOffset1,tOffset2,aMax,aMin)

sprintf('Fitting PVAJ model...')

%-- initialize global using parameters

% spatial parameters
u_azi = [0 45 90 135 180 225 270 315]';
u_ele = [-90 -45 0 45 90]';
s_data = [u_ele;u_azi]; % transform to this form for fitting
% time parameters
time = (1: nBins(1,1))' /nBins*(1500+tOffset1+tOffset2)/1000;
st_data = [u_ele;u_azi;time]; % transform to this form for fitting
% fitting initial parameters
mu = 0.91;
sig = sqrt(sqrt(2))/6;
baseline = spon;
acc_max = (aMax+tOffset1);
acc_min = (aMin+tOffset1);

% PSTH data
temporal_data = squeeze(mean(mean(PSTH_data(:,:,:),1),2)); % PSTH data according to time bin
spatial_data = permute(spatial_data,[2 1]);
y_data = permute(PSTH_data, [2 1 3]); % transform to azi*ele*timebin

i_gauss_time = pos_func([mu sig],time);
gauss_time = vel_func([mu sig],time);
d_gauss_time = acc_func([mu sig],time);
d2_gauss_time = jerk_func([mu sig],time);

u1 = i_gauss_time;
p1 = (u1'*gauss_time)/(u1'*u1);
u2 = gauss_time - p1*u1;
p21 = (u1'*d_gauss_time)/(u1'*u1);
p22 = (u2'*d_gauss_time)/(u2'*u2);
u3 = d_gauss_time - p21*u1 - p22*u2;
p31 = (u1'*d2_gauss_time)/(u1'*u1);
p32 = (u2'*d2_gauss_time)/(u2'*u2);
p33 = (u3'*d2_gauss_time)/(u3'*u3);
u4 = d2_gauss_time - p31*u1 - p32*u2 - p33*u3;

t_psth = y_data - baseline;
p_s_profile = zeros([length(u_azi), length(u_ele)]);
v_s_profile = zeros([length(u_azi), length(u_ele)]);
a_s_profile = zeros([length(u_azi), length(u_ele)]);
j_s_profile = zeros([length(u_azi), length(u_ele)]);

for j=1:length(u_ele),
    for i=1:length(u_azi),
        t_profile = squeeze(t_psth(i,j,:));
        coeff = (pinv([u1 u2 u3 u4])*squeeze(t_profile));
        p_s_profile(i,j) = coeff(1) - ...
            coeff(2)*p1 + ...
            coeff(3)*(-p21+p22*p1) + ...
            coeff(4)*(-p31+p32*p1+p33*(p21-p22*p1));
        
        v_s_profile(i,j) = coeff(2) - coeff(3)*p22 + ...
            coeff(4)*(-p32+p33*p22);
        
        a_s_profile(i,j) = coeff(3) - coeff(4)*p33;
        
        j_s_profile(i,j) = coeff(4);
    end
end

% normalise time and spatial profile
p_DC = (min(p_s_profile(:))+max(p_s_profile(:)))/2;
p_A = (max(p_s_profile(:))-min(p_s_profile(:)))/2;
p_space_profile = p_s_profile-p_DC;
p_space_profile = p_space_profile/p_A;

v_DC = (min(v_s_profile(:))+max(v_s_profile(:)))/2;
v_A = (max(v_s_profile(:))-min(v_s_profile(:)))/2;
v_space_profile = v_s_profile-v_DC;
v_space_profile = v_space_profile/v_A;

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

%% fit velocity spatial profile
%-- 1st, fit spatial profile1

LB = [0.001 0 -pi/2];
UB = [10 2*pi pi/2];

[~, max_idx] = max(v_space_profile(:));
[max_idx_a, max_idx_e] = ind2sub(size(v_space_profile), max_idx);

param = [0.01 u_azi(max_idx_a) u_ele(max_idx_e)];
recon_p = lsqcurvefit('cos_tuning', param,  s_data, ...
    p_space_profile(:), LB, UB, options);
p_n = recon_p(1);
p_a_0 = recon_p(2);
p_e_0 = recon_p(3);

param = [0.01 u_azi(max_idx_a) u_ele(max_idx_e)];
recon_v = lsqcurvefit('cos_tuning', param,  s_data, ...
    v_space_profile(:), LB, UB, options);
v_n = recon_v(1);
v_a_0 = recon_v(2);
v_e_0 = recon_v(3);

param = [0.01 u_azi(max_idx_a) u_ele(max_idx_e)];
recon_v = lsqcurvefit('cos_tuning', param,  s_data, ...
    a_space_profile(:), LB, UB, options);
a_n = recon_v(1);
a_a_0 = recon_v(2);
a_e_0 = recon_v(3);

param = [0.01 u_azi(max_idx_a) u_ele(max_idx_e)];
recon_v = lsqcurvefit('cos_tuning', param,  s_data, ...
    j_space_profile(:), LB, UB, options);
j_n = recon_v(1);
j_a_0 = recon_v(2);
j_e_0 = recon_v(3);


%%

p_DC = p_DC/v_A;
v_DC = v_DC/v_A;
a_DC = a_DC/a_A;
j_DC = j_DC/j_A;

%Fit linear parameters
A = p_A+v_A+a_A+j_A;
R_0 = baseline;
w_p = p_A/(p_A+v_A);
w_a = a_A/(p_A+v_A+a_A);
w_j = j_A/A;

%% fit PVAJ model
mu_t = mu;
sig_t = sig;

%Inital fits
param = [A, ...     %1
    R_0, ...   %2
    mu_t, ...  %3
    p_n, ...   %4
    p_a_0, ... %5
    p_e_0, ... %6
    p_DC, ...  %7
    v_n, ...   %8
    v_a_0, ... %9
    v_e_0, ... %10
    v_DC, ...  %11
    a_n, ...   %12
    a_a_0, ... %13
    a_e_0, ... %14
    a_DC, ...  %15
    j_n, ...   %16
    j_a_0, ... %17
    j_e_0, ... %18
    j_DC, ...  %19
    w_p, ...   %20
    w_a, ...   &21
    w_j];      %22

init_param = zeros(reps+1, length(param));
init_param(1,:) = param;

LB = [0.25*A, ...  %1  A
    0, ...           %2  R_0
    acc_max, ...           %3  mu_t
    0.001, ...       %4  p_n
    0, ...           %5  p_a_0
    -pi/2, ...       %6  p_e_0
    0, ...           %7  p_DC
    0.001, ...       %8  v_n
    0, ...           %9  v_a_0
    -pi/2, ...       %10 v_e_0
    0, ...           %11 v_DC
    0.001, ...       %12 a_n
    0, ...           %13 a_a_0
    -pi/2, ...       %14 a_e_0
    0, ...           %15 a_DC
    0.001, ...       %16 j_n
    0, ...           %17 j_a_0
    -pi/2, ...       %18 j_e_0
    0, ...           %19 j_DC
    0, ...           %20 w_p
    0, ...           %21 w_a
    0];              %22 w_j

UB = [4*A, ...     %1  A
    300, ...         %2  R_0
    acc_min, ...         %3  mu_t
    10, ...          %4  p_n
    2*pi, ...        %5  p_a_0
    pi/2, ...        %6  p_e_0
    1, ...           %7  p_DC
    10, ...          %8  v_n
    2*pi, ...        %9  v_a_0
    pi/2, ...        %10 v_e_0
    1, ...           %11 v_DC
    10, ...          %12 a_n
    2*pi, ...        %13 a_a_0
    pi/2, ...        %14 a_e_0
    1, ...           %15 a_DC
    10, ...          %16 j_n
    2*pi, ...        %17 j_a_0
    pi/2, ...        %18 j_e_0
    1, ...           %19 j_D
    1, ...           %20 w_p
    1, ...           %21 w_a
    1];              %22 w_j


rand_rss = zeros(reps+1,1);
rand_param = zeros(reps+1, length(param));
rand_jac = zeros(reps+1, length(param), length(param));

[rand_param(1,:),rand_rss(1),~,~,~,~,temp_jac] = lsqcurvefit('PVAJ_Model', ...
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
    
    [rand_param(ii,:),rand_rss(ii),~,~,~,~,temp_jac] = lsqcurvefit('PVAJ_Model', ...
        seed_param, st_data, y_data, LB, UB, options);
    rand_jac(ii,:,:) = full(temp_jac)'*full(temp_jac);
    
    if rand_rss(ii) < min_rss
        min_rss = rand_rss(ii);
        min_param = rand_param(ii,:);
    end
    
    
end

% find the best fit parameters according to rss
[~,min_inx] = min(rand_rss);
modelFitPara_PVAJ = rand_param(min_inx,:);
rss_PVAJ = rand_rss(min_inx);
jac_PVAJ = rand_jac(min_inx,:,:);

% calculate the final model fitting values
respon = PVAJ_Model(modelFitPara_PVAJ,st_data);
modelFitRespon_PVAJ = respon;

modelFit_PVAJ.V = VO_Model(modelFitPara_PVAJ(1:12),st_data);
modelFit_PVAJ.A = AO_Model(modelFitPara_PVAJ([1:4,13:20]),st_data);
modelFit_PVAJ.J = JO_Model(modelFitPara_PVAJ([1:4,21:28]),st_data);
modelFit_PVAJ.P = PO_Model(modelFitPara_PVAJ([1:4,29:end]),st_data);
%% analysis
data_num = 26*nBins;
para_num = 18;
BIC_PVAJ = BIC_fit(data_num,rss_PVAJ,para_num);
TSS = sum((PSTH_data(:) - mean(PSTH_data(:))).^2);
RSquared_PVAJ = 1 - rss_PVAJ/TSS;


end
