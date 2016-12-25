%-----------------------------------------------------------------------------------------------------------------------
%-- Direction2d_cue_conflict_fit.m -- Fits 2D cue conflict data with sum,
% product models, etc
%-- GCD and KM, starting 12/16/04
%-----------------------------------------------------------------------------------------------------------------------
function Direction2d_cue_conflict_fit(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

Path_Defs;
ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP

%get the column of values for azimuth and elevation and stim_type
temp_azimuth_moog = data.moog_params(HEADING,:,MOOG);
temp_azimuth_cam = data.moog_params(HEADING,:,CAMERAS);
preferred_azimuth = data.one_time_params(PREFERRED_AZIMUTH);
preferred_elevation = data.one_time_params(PREFERRED_ELEVATION);

%now, get the firing rates for all the trials 
temp_spike_rates = data.spike_rates(SpikeChan, :);                                                                                                                             

%get indices of any NULL conditions (for measuring spontaneous activity
null_trials = logical( (temp_azimuth_moog == data.one_time_params(NULL_VALUE)) & (temp_azimuth_cam == data.one_time_params(NULL_VALUE)) );

%now, remove trials from direction and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:length(temp_azimuth_moog);		% a vector of trial indices
select_trials= ( (trials >= BegTrial) & (trials <= EndTrial) );
azimuth_moog = temp_azimuth_moog(~null_trials & select_trials);
azimuth_cam = temp_azimuth_cam(~null_trials & select_trials);
% elevation = temp_elevation(~null_trials & select_trials);

spike_rates = temp_spike_rates(~null_trials & select_trials);

unique_azimuth_moog = munique(azimuth_moog');
unique_azimuth_cam = munique(azimuth_cam');
% unique_elevation = munique(elevation');


%% ADD CODE HERE FOR PLOTTING
%create basic matrix to represent each response vector

resp = [];
for i=1:length(unique_azimuth_moog)
    for j=1:length(unique_azimuth_cam)
            select = logical( (azimuth_moog==unique_azimuth_moog(i)) & (azimuth_cam==unique_azimuth_cam(j)) );
            if (sum(select) > 0)
               resp(j, i) = mean(spike_rates(select));
               resp_std(j,i) = std(spike_rates(select));
               resp_ste(j,i) =resp_std(j,i) / sqrt(length(find( (azimuth_moog==unique_azimuth_moog(i)) & (azimuth_cam==unique_azimuth_cam(j)) )) );
            end
     end
end 

global resp_cam_pad resp_ves_pad resp_conflict%need these global for the model fitting programs to run

% the first one is the conflict datase; the second turns the matrix of
% conflict response values into a long vector (for use in running the
% sequential F-test)
resp_conflict = resp(2:length(unique_azimuth_moog), 2:length(unique_azimuth_cam) );%8x8 matrix with visual axis on the rows (each row is a different visual stim value)
CONFLICT_STAT = [resp_conflict_norm(:,1);resp_conflict_norm(:,2);resp_conflict_norm(:,3);resp_conflict_norm(:,4);resp_conflict_norm(:,5);resp_conflict_norm(:,6);resp_conflict_norm(:,7);resp_conflict_norm(:,8)];

% the first is the camera (visual only - tuning) dataset; the second is the standard error
% of that dataset; the third is the camera dataset "padded" (8 identical
% vectors to create an 8x8 matrix)
resp_cam = resp( 2 : length(unique_azimuth_cam), 1 );%8x1 vector for visual information (column)
resp_cam_ste = resp_ste( 2 : length(unique_azimuth_cam), 1 );
resp_cam_pad = [resp_cam resp_cam resp_cam resp_cam resp_cam resp_cam resp_cam resp_cam];

% same idea for vestibular only dataset
resp_ves = resp( 1, 2:length(unique_azimuth_moog) );%1x8 vector with vestibular information (row)
resp_ves_ste = resp_ste( 1, 2:length(unique_azimuth_moog) );
resp_ves_pad = [resp_ves; resp_ves; resp_ves; resp_ves; resp_ves; resp_ves; resp_ves; resp_ves];


%------------------------------------------------------------------------
%ANOVAs
 
%perform a one-way ANOVA on visual response only (is the tuning
%significant?)
cam_select1 = logical ( (azimuth_moog==-9999) & (azimuth_cam==0) );
cam_select2 = logical ( (azimuth_moog==-9999) & (azimuth_cam==45) );
cam_select3 = logical ( (azimuth_moog==-9999) & (azimuth_cam==90) );
cam_select4 = logical ( (azimuth_moog==-9999) & (azimuth_cam==135) );
cam_select5 = logical ( (azimuth_moog==-9999) & (azimuth_cam==180) );
cam_select6 = logical ( (azimuth_moog==-9999) & (azimuth_cam==225) );
cam_select7 = logical ( (azimuth_moog==-9999) & (azimuth_cam==270) );
cam_select8 = logical ( (azimuth_moog==-9999) & (azimuth_cam==315) );

cam_spk1 = spike_rates(cam_select1);
cam_spk2 = spike_rates(cam_select2);
cam_spk3 = spike_rates(cam_select3);
cam_spk4 = spike_rates(cam_select4);
cam_spk5 = spike_rates(cam_select5);
cam_spk6 = spike_rates(cam_select6);
cam_spk7 = spike_rates(cam_select7);
cam_spk8 = spike_rates(cam_select8);

cam_stat = [cam_spk1' cam_spk2' cam_spk3' cam_spk4' cam_spk5' cam_spk6' cam_spk7' cam_spk8'];
anova1_cam = anova1(cam_stat,[],'off')

%perform a one-way ANOVA on vestibular response only
moog_select1 = logical ( (azimuth_moog==0) & (azimuth_cam==-9999) );
moog_select2 = logical ( (azimuth_moog==45) & (azimuth_cam==-9999) );
moog_select3 = logical ( (azimuth_moog==90) & (azimuth_cam==-9999) );
moog_select4 = logical ( (azimuth_moog==135) & (azimuth_cam==-9999) );
moog_select5 = logical ( (azimuth_moog==180) & (azimuth_cam==-9999) );
moog_select6 = logical ( (azimuth_moog==225) & (azimuth_cam==-9999) );
moog_select7 = logical ( (azimuth_moog==270) & (azimuth_cam==-9999) );
moog_select8 = logical ( (azimuth_moog==315) & (azimuth_cam==-9999) );

moog_spk1 = spike_rates(moog_select1);
moog_spk2 = spike_rates(moog_select2);
moog_spk3 = spike_rates(moog_select3);
moog_spk4 = spike_rates(moog_select4);
moog_spk5 = spike_rates(moog_select5);
moog_spk6 = spike_rates(moog_select6);
moog_spk7 = spike_rates(moog_select7);
moog_spk8 = spike_rates(moog_select8);

moog_stat = [moog_spk1' moog_spk2' moog_spk3' moog_spk4' moog_spk5' moog_spk6' moog_spk7' moog_spk8'];
anova1_moog = anova1(moog_stat,[],'off')

%perform a two-way ANOVA on visual + vestibular responses using all
%responses (does visual and/or vestibular stimulus value significantly
%contribute to variation in response?  is this purely additive, or is there
%significant interaction between vis and ves?)
moog_spk1_1 = spike_rates(logical ( (azimuth_cam==0) & (azimuth_moog==0) ));
moog_spk1_2 = spike_rates(logical ( (azimuth_cam==0) & (azimuth_moog==45) ));
moog_spk1_3 = spike_rates(logical ( (azimuth_cam==0) & (azimuth_moog==90) ));
moog_spk1_4 = spike_rates(logical ( (azimuth_cam==0) & (azimuth_moog==135) ));
moog_spk1_5 = spike_rates(logical ( (azimuth_cam==0) & (azimuth_moog==180) ));
moog_spk1_6 = spike_rates(logical ( (azimuth_cam==0) & (azimuth_moog==225) ));
moog_spk1_7 = spike_rates(logical ( (azimuth_cam==0) & (azimuth_moog==270) ));
moog_spk1_8 = spike_rates(logical ( (azimuth_cam==0) & (azimuth_moog==315) ));

moog_spk2_1 = spike_rates(logical ( (azimuth_cam==45) & (azimuth_moog==0) ));
moog_spk2_2 = spike_rates(logical ( (azimuth_cam==45) & (azimuth_moog==45) ));
moog_spk2_3 = spike_rates(logical ( (azimuth_cam==45) & (azimuth_moog==90) ));
moog_spk2_4 = spike_rates(logical ( (azimuth_cam==45) & (azimuth_moog==135) ));
moog_spk2_5 = spike_rates(logical ( (azimuth_cam==45) & (azimuth_moog==180) ));
moog_spk2_6 = spike_rates(logical ( (azimuth_cam==45) & (azimuth_moog==225) ));
moog_spk2_7 = spike_rates(logical ( (azimuth_cam==45) & (azimuth_moog==270) ));
moog_spk2_8 = spike_rates(logical ( (azimuth_cam==45) & (azimuth_moog==315) ));

moog_spk3_1 = spike_rates(logical ( (azimuth_cam==90) & (azimuth_moog==0) ));
moog_spk3_2 = spike_rates(logical ( (azimuth_cam==90) & (azimuth_moog==45) ));
moog_spk3_3 = spike_rates(logical ( (azimuth_cam==90) & (azimuth_moog==90) ));
moog_spk3_4 = spike_rates(logical ( (azimuth_cam==90) & (azimuth_moog==135) ));
moog_spk3_5 = spike_rates(logical ( (azimuth_cam==90) & (azimuth_moog==180) ));
moog_spk3_6 = spike_rates(logical ( (azimuth_cam==90) & (azimuth_moog==225) ));
moog_spk3_7 = spike_rates(logical ( (azimuth_cam==90) & (azimuth_moog==270) ));
moog_spk3_8 = spike_rates(logical ( (azimuth_cam==90) & (azimuth_moog==315) ));

moog_spk4_1 = spike_rates(logical ( (azimuth_cam==135) & (azimuth_moog==0) ));
moog_spk4_2 = spike_rates(logical ( (azimuth_cam==135) & (azimuth_moog==45) ));
moog_spk4_3 = spike_rates(logical ( (azimuth_cam==135) & (azimuth_moog==90) ));
moog_spk4_4 = spike_rates(logical ( (azimuth_cam==135) & (azimuth_moog==135) ));
moog_spk4_5 = spike_rates(logical ( (azimuth_cam==135) & (azimuth_moog==180) ));
moog_spk4_6 = spike_rates(logical ( (azimuth_cam==135) & (azimuth_moog==225) ));
moog_spk4_7 = spike_rates(logical ( (azimuth_cam==135) & (azimuth_moog==270) ));
moog_spk4_8 = spike_rates(logical ( (azimuth_cam==135) & (azimuth_moog==315) ));

moog_spk5_1 = spike_rates(logical ( (azimuth_cam==180) & (azimuth_moog==0) ));
moog_spk5_2 = spike_rates(logical ( (azimuth_cam==180) & (azimuth_moog==45) ));
moog_spk5_3 = spike_rates(logical ( (azimuth_cam==180) & (azimuth_moog==90) ));
moog_spk5_4 = spike_rates(logical ( (azimuth_cam==180) & (azimuth_moog==135) ));
moog_spk5_5 = spike_rates(logical ( (azimuth_cam==180) & (azimuth_moog==180) ));
moog_spk5_6 = spike_rates(logical ( (azimuth_cam==180) & (azimuth_moog==225) ));
moog_spk5_7 = spike_rates(logical ( (azimuth_cam==180) & (azimuth_moog==270) ));
moog_spk5_8 = spike_rates(logical ( (azimuth_cam==180) & (azimuth_moog==315) ));

moog_spk6_1 = spike_rates(logical ( (azimuth_cam==225) & (azimuth_moog==0) ));
moog_spk6_2 = spike_rates(logical ( (azimuth_cam==225) & (azimuth_moog==45) ));
moog_spk6_3 = spike_rates(logical ( (azimuth_cam==225) & (azimuth_moog==90) ));
moog_spk6_4 = spike_rates(logical ( (azimuth_cam==225) & (azimuth_moog==135) ));
moog_spk6_5 = spike_rates(logical ( (azimuth_cam==225) & (azimuth_moog==180) ));
moog_spk6_6 = spike_rates(logical ( (azimuth_cam==225) & (azimuth_moog==225) ));
moog_spk6_7 = spike_rates(logical ( (azimuth_cam==225) & (azimuth_moog==270) ));
moog_spk6_8 = spike_rates(logical ( (azimuth_cam==225) & (azimuth_moog==315) ));

moog_spk7_1 = spike_rates(logical ( (azimuth_cam==270) & (azimuth_moog==0) ));
moog_spk7_2 = spike_rates(logical ( (azimuth_cam==270) & (azimuth_moog==45) ));
moog_spk7_3 = spike_rates(logical ( (azimuth_cam==270) & (azimuth_moog==90) ));
moog_spk7_4 = spike_rates(logical ( (azimuth_cam==270) & (azimuth_moog==135) ));
moog_spk7_5 = spike_rates(logical ( (azimuth_cam==270) & (azimuth_moog==180) ));
moog_spk7_6 = spike_rates(logical ( (azimuth_cam==270) & (azimuth_moog==225) ));
moog_spk7_7 = spike_rates(logical ( (azimuth_cam==270) & (azimuth_moog==270) ));
moog_spk7_8 = spike_rates(logical ( (azimuth_cam==270) & (azimuth_moog==315) ));

moog_spk8_1 = spike_rates(logical ( (azimuth_cam==315) & (azimuth_moog==0) ));
moog_spk8_2 = spike_rates(logical ( (azimuth_cam==315) & (azimuth_moog==45) ));
moog_spk8_3 = spike_rates(logical ( (azimuth_cam==315) & (azimuth_moog==90) ));
moog_spk8_4 = spike_rates(logical ( (azimuth_cam==315) & (azimuth_moog==135) ));
moog_spk8_5 = spike_rates(logical ( (azimuth_cam==315) & (azimuth_moog==180) ));
moog_spk8_6 = spike_rates(logical ( (azimuth_cam==315) & (azimuth_moog==225) ));
moog_spk8_7 = spike_rates(logical ( (azimuth_cam==315) & (azimuth_moog==270) ));
moog_spk8_8 = spike_rates(logical ( (azimuth_cam==315) & (azimuth_moog==315) ));

conflict_stat1 = [moog_spk1_1' moog_spk1_2' moog_spk1_3' moog_spk1_4' moog_spk1_5' moog_spk1_6' moog_spk1_7' moog_spk1_8'];
conflict_stat2 = [moog_spk2_1' moog_spk2_2' moog_spk2_3' moog_spk2_4' moog_spk2_5' moog_spk2_6' moog_spk2_7' moog_spk2_8'];
conflict_stat3 = [moog_spk3_1' moog_spk3_2' moog_spk3_3' moog_spk3_4' moog_spk3_5' moog_spk3_6' moog_spk3_7' moog_spk3_8'];
conflict_stat4 = [moog_spk4_1' moog_spk4_2' moog_spk4_3' moog_spk4_4' moog_spk4_5' moog_spk4_6' moog_spk4_7' moog_spk4_8'];
conflict_stat5 = [moog_spk5_1' moog_spk5_2' moog_spk5_3' moog_spk5_4' moog_spk5_5' moog_spk5_6' moog_spk5_7' moog_spk5_8'];
conflict_stat6 = [moog_spk6_1' moog_spk6_2' moog_spk6_3' moog_spk6_4' moog_spk6_5' moog_spk6_6' moog_spk6_7' moog_spk6_8'];
conflict_stat7 = [moog_spk7_1' moog_spk7_2' moog_spk7_3' moog_spk7_4' moog_spk7_5' moog_spk7_6' moog_spk7_7' moog_spk7_8'];
conflict_stat8 = [moog_spk8_1' moog_spk8_2' moog_spk8_3' moog_spk8_4' moog_spk8_5' moog_spk8_6' moog_spk8_7' moog_spk8_8'];

conflict_stat = [conflict_stat1; conflict_stat2; conflict_stat3; conflict_stat4; conflict_stat5; conflict_stat6; conflict_stat7; conflict_stat8];
reps = length(moog_spk1_1);
anova2_all_conflict = anova2(conflict_stat, reps, 'hide')%rows are visual(cam), columns are vestibular(moog)

%-------------------------------------------------------------------------

% calculate mean spontaneous firing rate
spon_found = find(null_trials==1); 
spon_resp = mean(temp_spike_rates(spon_found));
spon_resp_pad = spon_resp*ones(8);

%--------------------------------------------------------------------------
%BRUTE FORCE METHOD OF MODEL FITTING - fit to mean responses without
%fmincon

% %find the value of Kp that leads to the least sum of square errors in a
% %weighted product model that uses mean responses to each stimulus condition
% for i=1:500,
%     squ_perror(i) = sum(sum((i/100 * resp_cam_pad.*resp_ves_pad + spon_resp_pad - resp_conflict).^2));
% end
% [min_perror,Ip]=min(squ_perror);
% Kp=Ip/100
% pmodel=Kp*resp_cam_pad.*resp_ves_pad + spon_resp_pad;
% 
% %figure for the weighted product model of best fit
% figure(2)
% Z1=pcolor(pmodel')
% shading interp
% title('Weighted Product Model of Best Fit')
% xlabel('Visual')
% ylabel('Vestibular')
% 
% %find the value of Kp that leads to the least mean square error in a
% %weighted product model that uses mean responses to each stimulus condition
% for i=1:100,
%    mean_squ_perror(i) = mean(mean((i/100 * resp_cam_pad.*resp_ves_pad + spon_resp_pad - resp_conflict).^2));
% end
% [mean_min_perror, Imp] = min(mean_squ_perror);
% Kmp=Imp/100;
% 
% %find the values of Ks and alpha that lead to the least sum of square
% %errors in a weighted sum model that uses mean responses to each stimulus
% %condition
% for i=1:100,
%     for j=1:100,
%         squ_serror(i,j)=sum(sum((i/100*(j/100*resp_cam_pad + (1-j/100)*resp_ves_pad) + spon_resp_pad - resp_conflict).^2));
%     end
% end
% [premin_serror,preIs]=min(squ_serror);
% [min_serror,Is]=min(premin_serror');
% alpha=Is/100
% Ks=(preIs(1,Is))/100
% smodel=Ks * (alpha*resp_cam_pad + (1-alpha)*resp_ves_pad) + spon_resp_pad;
% 
% figure(3)
% Z2=pcolor(smodel')
% shading interp
% title('Weighted Sum Model of Best Fit')
% xlabel('Visual')
% ylabel('Vestibular')
% 
% %find chi2 values for both models; if n=63, critical chi2 = 82.53 --> if
% %chi2 > 82.53, we can reject the null hypothesis (the two are significantly
% %different)
% %expected = pmodel or smodel
% %observed = resp_conflict
% chi2product = sum(sum(((resp_conflict-pmodel).^2)./pmodel))
% chi2sum = sum(sum(((resp_conflict-smodel).^2)./smodel))
% 
% %find the values of Ks and alpha that lead to the least mean square error
% %in a weighted sum model that uses mean responses to each stimulus
% %condition
% for i=1:100,
%     for j=1:100,
%         mean_squ_serror(i,j)=mean(mean((i/100*(j/100*resp_cam_pad + (1-j/100)*resp_ves_pad) + spon_resp_pad - resp_conflict).^2));
%     end
% end
% [mean_premin_serror, preIms]=min(mean_squ_serror);
% [mean_min_serror, Ims]=min(mean_premin_serror');
% alpham=Ims/100;
% Kms=(preIms(1,Ims))/100;

%------------------------------------------------------------------
%FMINCON METHOD FOR MODEL FITTING

%use fmincon to fit models to the data -- mean data

%generate an initial guess
q = [1 1];

%set up parameters for fmincon
A=[]; B=[]; Aeq=[]; Beq=[];
LB = [-50 0]; UB = [50 5]; NONLCON=[];

OPTIONS = OPTIMSET('fmincon');
OPTIONS = OPTIMSET('Largescale', 'off', 'LevenbergMarquardt', 'on', 'MaxFunEvals', 5000, 'MaxIter', 5000);

%run fmincon
N_reps = 10;
wiggle = 0.8;
testpars = []; err=[];
for j=1:N_reps
    rand_factor = rand(2,1) * wiggle + (1-wiggle/2);
    temp_q = q'.* rand_factor;
    testpars{j} = fmincon('prodmoderr',temp_q,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS);
    err(j) = prodmoderr(testpars{j});%all of the values in these keep coming back equal
end

%now that we have models and errors, find the best fit
[min_err min_err_indx] = min(err)
pars = testpars{min_err_indx}

%create a matrix containing the weighted product model
PRODMOD = prodmod(resp_cam_pad,resp_ves_pad,pars);

%create a column vector, in which the columns from PRODMOD are stacked on
%each other into a single column (to be compared to CONFLICT_STAT)
PRODMOD_STAT = [PRODMOD(:,1);PRODMOD(:,2);PRODMOD(:,3);PRODMOD(:,4);PRODMOD(:,5);PRODMOD(:,6);PRODMOD(:,7);PRODMOD(:,8)];

%perform chi square analysis of model versus actual mean response data
%by comparing expected (PRODMOD_STAT) to actual (CONFLICT_STAT)
%PRODMOD_CHI2 = sum(sum(((CONFLICT_STAT - PRODMOD_STAT).^2)/PRODMOD_STAT))
%with 64 values to compare, n=63, and the critical chi square value is
%82.53 --> if chi square > 82.53, we can reject the null hypothesis and
%find that the model and the actual data are significantly different

%calculate some statistics for the regression between the model and the actual mean response data, then
%plot the scatter of model versus actual mean response data
% PSTATS = regstats(PRODMOD_STAT,CONFLICT_STAT,'linear',{'R','beta','covb','yhat','r','mse'})
% P_R2 = (PSTATS.R)^2
% P_beta = PSTATS.beta
% P_cov = PSTATS.covb
% figure(3)
% scatter(PRODMOD_STAT,CONFLICT_STAT)
% title('Regression Analysis of Product Model Fit')
% xlabel('Product Model Values')
% ylabel('Actual Values')


%--------------------------------------------------------------------------

%use fmincon to generate a weighted sum model

%generate an initial guess
r = [0.1 1 1];%can't set an entry to zero, or it will be zero everytime in the fmincon loop

%set up parameters for fmincon
A=[]; B=[]; Aeq=[]; Beq=[];
LB = [-50 0 0]; UB = [50 5 1]; NONLCON=[];

OPTIONS = OPTIMSET('fmincon');
OPTIONS = OPTIMSET('Largescale', 'off', 'LevenbergMarquardt', 'on', 'MaxFunEvals', 5000, 'MaxIter', 5000);

%run fmincon
N_reps = 10;
wiggle = 0.5;
s_testpars = []; s_err=[];
for j=1:N_reps
    rand_factor = rand(3,1) * wiggle + (1-wiggle/2);
    temp_r = r'.* rand_factor;
    s_testpars{j} = fmincon('summoderr',temp_r,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS);
    s_err(j) = summoderr(s_testpars{j});
end

%now that we have models and errors, find the best fit
[min_s_err min_s_err_indx] = min(s_err)
s_pars = s_testpars{min_s_err_indx}

%create a matrix containing the weighted sum model
SUMMOD = summod(resp_cam_pad,resp_ves_pad,s_pars);

%create a column vector for this model similar to PRODMOD_STAT
SUMMOD_STAT = [SUMMOD(:,1);SUMMOD(:,2);SUMMOD(:,3);SUMMOD(:,4);SUMMOD(:,5);SUMMOD(:,6);SUMMOD(:,7);SUMMOD(:,8)];

%calculate chi square value comparing the model (SUMMOD_STAT) to the actual
%data (CONFLICT_STAT)
%SUMMOD_CHI2 = sum(sum(((CONFLICT_STAT - SUMMOD_STAT).^2)/SUMMOD_STAT))

%calculate some statistics for the regression between the model and the actual mean response data, then
%plot the scatter of model versus actual mean response data
% SSTATS = regstats(SUMMOD_STAT,CONFLICT_STAT,'linear',{'R','beta','covb','yhat','r','mse'})
% S_R2 = (SSTATS.R)^2
% S_beta = SSTATS.beta
% S_cov = SSTATS.covb
% figure(5)
% scatter(SUMMOD_STAT,CONFLICT_STAT)
% title('Regression Analysis of Sum Model Fit')
% xlabel('Sum Model Values')
% ylabel('Actual Values')

%--------------------------------------------------------------------------
%plot the best fit models

figure(2)
Z1=pcolor((PRODMOD)');
%prodmod produces a matrix like resp_conflict, where row indicates visual
%stimulus value and column indicates vestibular stimulus value --> need to flip it to have the visual on the horizontal 
shading interp
title('Weighted Product Model')
xlabel('Visual')
ylabel('Vestibular')


figure(3)
Z2=pcolor((SUMMOD)');
%summod produces a matrix like resp_conflict, where row indicates visual
%stimulus value and column indicates vestibular stimulus value --> need to flip it to have the visual on the horizontal 
shading interp
title('Weighted Sum Model')
xlabel('Visual')
ylabel('Vestibular')

%-------------------------------------------------------------------------

%use fmincon to generate a combination model (product + sum elements)

%generate an initial guess
K = [1 1 0.5 0.1];%can't set an entry to zero, or it will be zero everytime in the fmincon loop

%set up parameters for fmincon
A=[]; B=[]; Aeq=[]; Beq=[];
LB=[0 0 0 -50]; UB=[5 5 1 50]; NONLCON=[];

OPTIONS = OPTIMSET('fmincon');
OPTIONS = OPTIMSET('Largescale', 'off', 'LevenbergMarquardt', 'on', 'MaxFunEvals', 5000, 'MaxIter', 5000);

%run fmincon
N_reps = 10;
wiggle = 0.5;
c_testpars = []; c_err=[];
for j=1:N_reps
    rand_factor = rand(4,1) * wiggle + (1-wiggle/2);
    temp_K = K'.* rand_factor;
    c_testpars{j} = fmincon('combomoderr',temp_K,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS);
    c_err(j) = combomoderr(c_testpars{j});
end

%now that we have models and errors, find the best fit
[min_c_err min_c_err_indx] = min(c_err)
c_pars = c_testpars{min_c_err_indx}

%create a matrix containing the weighted sum model - normalized
%SUMMOD = summod(resp_cam_pad_norm,resp_ves_pad_norm,s_pars);

%create a matrix containing the weighted sum model - pspont
COMBOMOD = combomod(resp_cam_pad,resp_ves_pad,c_pars);

%create a column vector for this model similar to PRODMOD_STAT
COMBOMOD_STAT = [COMBOMOD(:,1);COMBOMOD(:,2);COMBOMOD(:,3);COMBOMOD(:,4);COMBOMOD(:,5);COMBOMOD(:,6);COMBOMOD(:,7);COMBOMOD(:,8)];

%plot the combo model
figure (4)
Z3=pcolor((COMBOMOD)');
shading interp
title ('Combination Model')
xlabel ('Visual')
ylabel ('Vestibular')

%--------------------------------------------------------------------------
%SEQUENTIAL F-TEST COMPARING MODEL FITS TO EACH OTHER

plow = 2;
pmid = 3;
phigh = 4;
nn = 64;
sslow = min_err;
ssmid = min_s_err;
sshigh = min_c_err;

Fseq_sp = ((sslow - ssmid)/(pmid - plow))/(ssmid/(nn - pmid));
Fp_sp = 1 - fcdf(Fseq_sp,(pmid - plow), (nn - pmid))

Fseq_pc = ((sslow - sshigh)/(phigh - plow))/(sshigh/(nn - phigh));
Fp_pc = 1 - fcdf(Fseq_pc,(phigh - plow), (nn - phigh))

Fseq_sc = ((ssmid - sshigh)/(phigh - pmid))/(sshigh/(nn - phigh));
Fp_sc = 1 - fcdf(Fseq_sc,(phigh - pmid), (nn - phigh))

