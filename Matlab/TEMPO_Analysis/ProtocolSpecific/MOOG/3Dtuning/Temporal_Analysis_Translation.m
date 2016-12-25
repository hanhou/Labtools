%----------------------------------------------------------------------------------------------------------------------
%Frequency_Analysis.m calculates the DFTR, phase, latency, velocity and acceleration
%components, preferred directions (for MFR, velocity and acceleration),
%DDI's, vector sum amplitudes, correlation of 3d tuning and anova for neuronal responses IN THE
%VESTIBULAR (stim_type = 1) CONDITION.

function Frequency_Analysis(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
Path_Defs;
ProtocolDefs;  %contains protocol specific keywords - 1/4/01 BJP

%get the column of values for azimuth and elevation and stim_type
temp_azimuth = data.moog_params(AZIMUTH,:,MOOG);
temp_elevation = data.moog_params(ELEVATION,:,MOOG);
temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG);
temp_amplitude = data.moog_params(AMPLITUDE,:,MOOG);
% fixation data added by CRF, 3/2007
temp_fix_x    =  data.moog_params(FIX_X,:,MOOG);
temp_fix_y    =  data.moog_params(FIX_Y,:,MOOG);
temp_fix_x(isnan(temp_fix_x)) = 0;
temp_fix_y(isnan(temp_fix_y)) = 0;

temp_spike_data = data.spike_data(1,:);
dummy_spike_data = temp_spike_data;
%mean firing rate of each trial depending on the start and stop offsets
temp_spike_rates = data.spike_rates(SpikeChan, :);

%get indices of any NULL conditions (for measuring spontaneous activity
null_trials = logical( (temp_azimuth == data.one_time_params(NULL_VALUE)) );

%now, remove trials from direction and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:length(temp_azimuth); % a vector of trial indices
%at any given point, cell response cannot be greater than 1
abnormal = find(temp_spike_data > 1);
temp_spike_data(1,abnormal) = 1;
bad_trials = find(temp_spike_rates > 3000);   % cut off 3k frequency which definately is not cell's firing response

if ( bad_trials ~= NaN)
    select_trials= ( (trials >= BegTrial) & (trials <= EndTrial) & (trials~=bad_trials) );
else 
    select_trials= ( (trials >= BegTrial) & (trials <= EndTrial) ); 
end

azimuth = temp_azimuth(~null_trials & select_trials);
elevation = temp_elevation(~null_trials & select_trials);
stim_type = temp_stim_type(~null_trials & select_trials);
fix_x = temp_fix_x(~null_trials & select_trials);
fix_y = temp_fix_y(~null_trials & select_trials);
spike_rates = temp_spike_rates(~null_trials & select_trials);

% **** Not used anywhere else in the code; commenting this out **** -CRF 2/2007
% %baseline firing rate
% spon_found = find(null_trials==1);
% spon_resp = mean(temp_spike_rates(null_trials));
% for i = 1:length(spon_found)
%     total_null(i) = temp_spike_rates(spon_found(i));
% end
% spon_count = floor(sum(temp_spike_rates(null_trials))/20);
%
% total_trials = floor( length(azimuth) /(26*length(unique_condition_num)) );

unique_azimuth = munique(azimuth');
unique_elevation = munique(elevation');
unique_stim_type = munique(stim_type');
unique_fix_x    =  munique(fix_x');
unique_fix_y    =  munique(fix_y');    

temp_condition_num = temp_stim_type;
condition_num = stim_type;
h_title{1}='Vestibular';
h_title{2}='Visual';
h_title{3}='Combined';
unique_condition_num = munique(condition_num');

timebin=50;
% sample frequency depends on test duration
frequency=length(temp_spike_data)/length(select_trials); % actually its 5000 now for 5 sec 
% length of x-axis
x_length = frequency/timebin;
% x-axis for plot PSTH
x_time=1:(frequency/timebin);

% find spontaneous trials which azimuth,elevation,stim_type= 9999
stim_duration = length(temp_spike_data)/length(temp_azimuth);

% remove null trials, bad trials, and trials outside Begtrial~Endtrial
Discard_trials = find(null_trials==1 | trials <BegTrial | trials >EndTrial);
for i = 1 : length(Discard_trials)
    temp_spike_data( 1, ((Discard_trials(i)-1)*stim_duration +1) :  Discard_trials(i)*stim_duration ) = 9999;
end
spike_data = temp_spike_data( temp_spike_data~=9999 );
% count spikes from raster data (spike_data)
max_count = 1;
time_step = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
current_stim = 1;  % vestibular=1, visual= 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k=1:length(unique_condition_num)  % it's the stim_type
    for j=1:length(unique_elevation)
        for i=1: length(unique_azimuth)         % SELECT now restricted to only 0 deg fixation trials -- CRF 3/2007
            select = logical(azimuth==unique_azimuth(i) & elevation==unique_elevation(j) & condition_num==unique_condition_num(k) & fix_x==0 & fix_y==0 );            
            if sum(select) > 0
                resp{k}(j,i) = mean(spike_rates(select));
                act_found = find( select==1 );
                % count spikes per timebin on every same condition trials
                clear temp_count dummy_count;
                for repeat=1:length(act_found) 
                    for n=1:(x_length)
                        temp_count(repeat,n)=sum(spike_data(1,(frequency*(act_found(repeat)-1)+time_step):(frequency*(act_found(repeat)-1)+n*timebin)));
                        dummy_count{repeat}(n) = sum(spike_data(1,(frequency*(act_found(repeat)-1)+time_step):(frequency*(act_found(repeat)-1)+n*timebin)));
                        time_step=time_step+timebin;
                    end
                    time_step=1;                    
                    
                    if k == current_stim
                        count_condition{i,j,repeat} = dummy_count{repeat};
                    end
                end
                
                count_y_trial{i,j,k}(:,:) = temp_count;  % each trial's PSTH in vestibular condition
                if k == 1
                    count_y_trial1{i,j}(:,:) = temp_count;
                end
                
                dim=size(temp_count);
                if dim(1) > 1;
                    count_y{i,j,k} = mean(temp_count);
                else
                    count_y{i,j,k}= temp_count;     % for only one repetition cases
                end
            else
                resp{k}(j,i) = 0; 
                count_y{i,j,k}=0;
            end   
            % normalize count_y
            if max(count_y{i,j,k})~=0;
                count_y_norm{i,j,k}=count_y{i,j,k} / max(count_y{i,j,k});
            else
                count_y_norm{i,j,k}=0;
            end
        end
    end
    % now find the peak
    [row_max, col_max] = find( resp{k}(:,:)==max(max(resp{k}(:,:))) );
    % it is likely there are two peaks with same magnitude, choose the first one arbitraly
    row_m{k}=row_max(1);
    col_m{k}=col_max(1);
    if max(count_y{col_max(1), row_max(1), k})~=0;
        count_y_max{k} = count_y{col_max(1), row_max(1), k} / max(count_y{col_max(1), row_max(1), k});
    else
        count_y_max{k} =0;
    end
    % find the largest y to set scale later
    if max(count_y{col_max(1), row_max(1), k}) > max_count
        max_count = max(count_y{col_max(1), row_max(1), k});
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
row_m{1};
col_m{1};

% order in which the directions are plotted
plot_col = [1 1 1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 6 6 6 7 7 7 8 8 8];
plot_row = [5 4 3 2 1 4 3 2 4 3 2 4 3 2 4 3 2 4 3 2 4 3 2 4 3 2];

% find empty trials
sum_dat = zeros(1,100);
for j = 1:26
    sum_dat = sum_dat + count_y{plot_col(j), plot_row(j), current_stim};
end
% set length of data to read from PSTH
for ii = 1:100
    if sum(sum_dat(ii:end)) == 0
        span = ii;
        break
    end
end
if span < 80
    span = 68; % usually this is how much data there is
else % this is one of katsu's cells
    span = 80;
end

% perform frequency analysis
for j=1:26
    gdat = count_y{plot_col(j), plot_row(j), current_stim};
    temp_dat{j} = gdat;
    temp_dat1{j} = gdat;
    %---------------------------------------
    % find hilbert transform to get peak.
    shift_win(j) = 0;
    peak_time(j) = 40;
    flagged_dir(j) = 0;
    % running mean
    binsize = 3; % take mean of 3 neighbour bins
    slidesize = 1; % with a slide of 1 bin
    nn = 1;
    start_bin = 1;
    end_bin = start_bin + binsize -1;
    while end_bin < 100 % the length of gdat
        if sum( gdat(start_bin : end_bin) )~=0 % avoid divided by 0 later          
            runmean(nn) = mean( gdat(start_bin : end_bin) );
        else 
            runmean(nn) = 0;
        end
        start_bin = start_bin + slidesize;
        end_bin = start_bin + binsize -1;
        nn = nn + 1;
    end
    hildata = zeros(1,100);
    hildata(2:length(runmean)+1) = runmean;
    dummy_dat{j} = hildata(18:63) - mean(hildata(14:20));  % WHY SUBRACT THE MEAN OF THIS (ARBITRARY?) 350 MS THAT OVERLAPS WITH THE 1.5 SEC ANALYSIS WINDOW? -CRF
    total_data{j} = gdat(1:span);        
    hil{j} = hilbert(dummy_dat{j});% find the peak within the window that we are considering (1.5 sec)
    abb = abs(hil{j});
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% added by
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Katsu
    % Here, pull out Hilbert {J}, but there is 46 points
    hhh{j} = abb;
    Hilbert{j}=hhh{j}(4:end-3);%reduce to 40   by Katsu
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    HT = abs(hil{j}) - mean(abb(4:5));% remove the mean to analyze just the peak or else the center of mass will be biased
    HT(HT<0) = 0;
    %find four largest bins to find the peak of the envelope. The bin
    %with the largest mean over 3 bins will be the peak.
    dum_HT = HT(13:37);
    mm = find(dum_HT == max(dum_HT));
    m(1) = mm(1);
    dum_HT(m(1)) = 0;
    mm = find(dum_HT == max(dum_HT));
    m(2) = mm(1);
    dum_HT(m(2)) = 0;
    mm = find(dum_HT == max(dum_HT));
    m(3) = mm(1);
    dum_HT(m(3)) = 0;
    mm = find(dum_HT == max(dum_HT));
    m(4) = mm(1);
    dum_HT(m(4)) = 0;
    dum_HT = HT(13:37);
    for i = 1:4
        if m(i) == length(dum_HT)
            me(i) = mean(dum_HT(m(i)-1:m(i)));
        elseif m(i) == 1;
            me(i) = mean(dum_HT(m(i):m(i)+1));
        else
            me(i) = mean(dum_HT(m(i)-1:m(i)+1));
        end
    end
    maxbin = find(me == max(me));
    peak_temp = m(maxbin)+29;
    
    HT = HT(4:end-3);
    H_time = 1:.05:3;
    for hh = 1:length(HT)
        CMsum(hh) = HT(hh)*H_time(hh);
    end
    if sum(HT) > 0
        center_mass = sum(CMsum)/sum(HT);
        peak_t = find( abs(H_time - center_mass) == min(abs(H_time - center_mass)) );% point which is closest to the center of mass is used.
        
        if (peak_t+20)*.05 >= 1.5 & (peak_t+20)*.05 <= 2.7
            peak_time(j) = peak_t(1) + 20; %if the center of mass lies right in the middle of two points, arbitrarily select the first point.
        else 
            peak_time(j) = peak_temp(1);
            flagged_dir(j) = 1;
        end
    end
    
    % stim peak time = 2 sec = 40th bin     % DOESN'T THIS ASSUME A VELOCITY RESPONSE? -CRF
    shift_win(j) = peak_time(j) - 40; %shift window for fourier transform to get rid of phase shift caused by delayed response.
    span_dum(j) = span;       
    %if the window is being shifted to extend into the region where the
    %recording has stopped, then we will randomly sample from the front of
    %the response and fill out the window.
    signal = total_data{j}(14:20);
    ran_in{j} = randperm(length(signal));
    rand_sig = signal(ran_in{j});
    span_new(j) = span;
    
    if shift_win(j) > 8 & span == 68 % if we shift by more than 8 bins, there is not recorded data at the end of the window
        for i = 1:(shift_win(j) - 8)
            total_data{j}(span+i) = rand_sig(i);
            span_new(j) = span+i;
        end
    end
    
    for k=1:40
        gaussdata(k)=total_data{j}(k+20+shift_win(j)); %window is shifted depending on how delayed the peak of the envelope (from the hilbert trans) is
    end
    gauss_dat{j} = gaussdata;
    gaustime=[x_time(20+shift_win(j)):x_time(59+shift_win(j))];
    gausstime=0.05*gaustime;    
    %calculate DFT ratio
    %40 pt DFT with mean removed
    [f, amp, resp_phase] = FT(gausstime, gauss_dat{j}, 40, 1, 0);
    f1 = mean(amp(2:4));
    f2 = mean(amp(5:end));
    if f2 == 0
        f_rat = 0;
    else
        %DFTR
        f_rat = f1/f2;
    end
    %phase of response
    max_p = mean(resp_phase(find(amp == max(amp(2:4)))));
    max_phase(j) = max_p(1);
    %calculate p-values
    [dft_p(j) dft_p_vel(j) dft_p_acc(j)] = DFT_perm(gauss_dat{j}, f_rat, max_phase(j), 1000);
    fourier_ratio(j) = f_rat;
    mean_rates(j) = sum(temp_dat{j}(30:50));
end

temp_max = find(fourier_ratio == max(fourier_ratio));
shift_phase = max_phase; % shift negative angle by 360
angle_add = zeros(1,length(shift_phase));
angle_add(find(shift_phase<0)) = 2*pi;
shift_phase = shift_phase+angle_add;
m_phase = shift_phase(temp_max);

% calculate the preferred direction using vector sum
DFT_ratio = zeros(5,8);
for i = 1:26
    DFT_ratio(plot_row(i), plot_col(i)) = fourier_ratio(i);
end

%----------------------------------------------------------------------------
%calculate preferred direction for velocity and acceleration by resolving
%the phase into two components - 180/0 and -90/90
vel_sum = zeros(5,8);
acc_sum = zeros(5,8);
for i = 1:26
    pol = -1;%since we want polarity to be positive for phases near 180 and negative for phases near 0
    temp_phase = max_phase(i);
    fourier_ratio_vel(i) = pol*fourier_ratio(i)*cos(temp_phase);
    fourier_ratio_acc(i) = fourier_ratio(i)*sin(temp_phase);
    vel_sum(plot_row(i), plot_col(i)) = pol*fourier_ratio(i)*cos(temp_phase);
    acc_sum(plot_row(i), plot_col(i)) = fourier_ratio(i)*sin(temp_phase);
end

[mfr_pref_az mfr_pref_el mfr_pref_amp] = vectorsum(resp{current_stim});
[pref_az_full pref_el_full pref_amp_full] = vectorsum(DFT_ratio); % ADDED pref dir based on full (raw) DFTR  -CRF
[pref_az_vel pref_el_vel pref_amp_vel] = vectorsum(vel_sum);
[pref_az_acc pref_el_acc pref_amp_acc] = vectorsum(acc_sum);
%----------------------------------------------------------------------------

%find the dft ratio of each individual trial
trial_time = 0:.05:5;
Azims = [0 45 90 135 180 225 270 315];
Elevs = [-90 -45 0 45 90];
azimuth1 = azimuth(find(condition_num == 1 & fix_x == 0 & fix_y == 0));  % restricted to only 0 deg fixation trials -- CRF 3/2007
elevation1 = elevation(find(condition_num == 1 & fix_x == 0 & fix_y == 0));
spike_rates_dum = spike_rates(find(condition_num == 1 & fix_x == 0 & fix_y == 0));
temp_spike_data = data.spike_data(1,:);

% puts the binned spike count data for each trial in a cell array (all_count)
ii = 1;
for i = 1:5000:length(temp_spike_data)  
   tdat = temp_spike_data(i:i+4999);
   step = 1:50:5001;
   for n=1:(x_length)
       all_count{ii}(n) = sum(tdat(step(n):step(n+1)-1));
   end
   ii = ii+1;
end

% For SELECT, changed azimuth/elevation, etc. to temp_azimuth/temp_elevation to include null trials.
% This was necessary to match the trial numbers for the parsing-by-repetition in the ANOVA code below. -CRF 3/2007
for k=1:length(unique_condition_num)
    for j=1:length(unique_elevation)
        for i=1: length(unique_azimuth)          % SELECT now restricted to only 0 deg fixation trials -CRF 3/2007
            select = logical( (temp_azimuth==unique_azimuth(i)) & (temp_elevation==unique_elevation(j)) & (temp_condition_num==unique_condition_num(k)) & (temp_fix_x==0) & (temp_fix_y==0) );
            % get rid off -90 and 90 cases
            if (sum(select) > 0)
                act_found = find( select==1 );
                % the following 4 lines are the same for each trial, so I moved them outside the FOR loop below  -CRF 3/2007
                current_azi_ele = find( (Azims(plot_col) == unique_azimuth(i)) & (Elevs(plot_row) == unique_elevation(j)) );
                trial_peak = peak_time(current_azi_ele);
                shift_trial = shift_win(current_azi_ele);
                trial_span = span;
                for repeat=1:length(act_found) 
                    temp_trial_data = all_count{act_found(repeat)};%changes for each trial
                    signal = temp_trial_data(14:20);%changes for each trial
                    rand_sig = signal(ran_in{current_azi_ele});%changes for each trial
                    temp_trial_data = temp_trial_data(1:trial_span);%changes for each trial
                    
                    if shift_trial > 8 & span == 68 % if we shift by more than 8 bins, there is not recorded data at the end of the window
                        for ii = 1:(shift_trial - 8)
                            temp_trial_data(trial_span+ii) = rand_sig(ii);
                        end
                    end
                    % now we have to calculate the dft ratio for each trial
                    [trial_f trial_amp tri_phase] = FT(trial_time((21+shift_trial):(60+shift_trial)), temp_trial_data((21+shift_trial):(60+shift_trial)), 40, 1, 0);
                    f1_trial = mean(trial_amp(2:4));
                    f2_trial = mean(trial_amp(5:end));
                    
                    if f2_trial == 0
                        dft_trial_all(act_found(repeat)) = 0;
                    else
                        dft_trial_all(act_found(repeat)) = f1_trial/f2_trial;
                    end
                    trial_phase_all(act_found(repeat)) = mean(tri_phase(find(trial_amp == max(trial_amp(2:4)))));
                    dft_velocity_all(act_found(repeat)) = -dft_trial_all(act_found(repeat))*cos(trial_phase_all(act_found(repeat)));
                    dft_acceleration_all(act_found(repeat)) = dft_trial_all(act_found(repeat))*sin(trial_phase_all(act_found(repeat)));
                end
            end
        end
    end
end

% repeats all of the above, but for only VESTIBULAR trials (via azimuth1 and elevation1)
% (only other difference I can find is referenced by '***' below  -CRF)
for j=1:length(unique_elevation)
    for i=1: length(unique_azimuth)
        select = logical( (azimuth1==unique_azimuth(i)) & (elevation1==unique_elevation(j)) );
        % get rid off -90 and 90 cases
        if (sum(select) > 0)
            act_found = find( select==1 );
            % the following 4 lines are the same for each trial, so I moved them outside the FOR loop below  -CRF 3/2007
            current_azi_ele = find( (Azims(plot_col) == unique_azimuth(i)) & (Elevs(plot_row) == unique_elevation(j)) );
            trial_peak = peak_time(current_azi_ele);                
            shift_trial = shift_win(current_azi_ele);            
            trial_span = span;                      
            for repeat=1:length(act_found)
                temp_trial_data = count_y_trial1{i,j}(repeat,:);%changes for each trial   ***uses count_y_trial1 instead of all_count, but the result is the same!  -CRF
                signal = temp_trial_data(14:20);%changes for each trial
                rand_sig = signal(ran_in{current_azi_ele});%changes for each trial
                temp_trial_data = temp_trial_data(1:trial_span);%changes for each trial
                
                if shift_trial > 8 & span == 68 % if we shift by more than 8 bins, there is not recorded data at the end of the window
                    for ii = 1:(shift_trial - 8)
                        temp_trial_data(trial_span+ii) = rand_sig(ii);
                    end
                end
                % now we have to calculate the dft ratio for each trial
                [trial_f trial_amp tri_phase] = FT(trial_time((21+shift_trial):(60+shift_trial)), temp_trial_data((21+shift_trial):(60+shift_trial)), 40, 1, 0);
                f1_trial = mean(trial_amp(2:4));
                f2_trial = mean(trial_amp(5:end));

                if f2_trial == 0
                    dft_trial(act_found(repeat)) = 0;
                else
                    dft_trial(act_found(repeat)) = f1_trial/f2_trial;
                end
                trial_phase(act_found(repeat)) = mean(tri_phase(find(trial_amp == max(trial_amp(2:4)))));
                dft_velocity(act_found(repeat)) = -dft_trial(act_found(repeat))*cos(trial_phase(act_found(repeat)));
                dft_acceleration(act_found(repeat)) = dft_trial(act_found(repeat))*sin(trial_phase(act_found(repeat))); 
            end
        end
    end
end


% %---------------------------------------------------------------------------------------
% % The following is the old ANOVA code.  Commenting out in favor of Narayan's version below
% %---------------------------------------------------------------------------------------
%  
% % Fixed to be correct for Vary_Fixation blocks -CRF 3/2007
% if length(unique_fix_x) >= length(unique_fix_y)
%     num_eye_pos = length(unique_fix_x);
% else
%     num_eye_pos = length(unique_fix_y);
% end
% 
% trials_per_rep = 26 * length(unique_condition_num) * length(unique_fix_x) * length(unique_fix_y) + num_eye_pos; % the last term corresponds to the null trial(s), one for each eye position
% repetitions = floor( (EndTrial-(BegTrial-1)) / trials_per_rep);
% 
% % also need to correct the length of the dft vectors, which were truncated at the last relevant-condition trial
% dft_trial_all(end+1 : repetitions*trials_per_rep) = 0;
% trial_phase_all(end+1 : repetitions*trials_per_rep) = 0;
% dft_velocity_all(end+1 : repetitions*trials_per_rep) = 0;
% dft_acceleration_all(end+1 : repetitions*trials_per_rep) = 0;
% 
% % first parse raw data into repetitions, ***including null trials***
% for q = 1:repetitions
%     azimuth_rep{q} = temp_azimuth(trials_per_rep*(q-1)+BegTrial : trials_per_rep*q+BegTrial-1);
%     elevation_rep{q} = temp_elevation(trials_per_rep*(q-1)+BegTrial : trials_per_rep*q+BegTrial-1);
%     condition_num_rep{q} = temp_stim_type(trials_per_rep*(q-1)+BegTrial : trials_per_rep*q+BegTrial-1);
%     fix_x_rep{q} = temp_fix_x(trials_per_rep*(q-1)+BegTrial : trials_per_rep*q+BegTrial-1);      % Added by CRF 3/2007 
%     fix_y_rep{q} = temp_fix_y(trials_per_rep*(q-1)+BegTrial : trials_per_rep*q+BegTrial-1);      % Added by CRF 3/2007
%     spike_rates_rep{q} = temp_spike_rates(trials_per_rep*(q-1)+BegTrial : trials_per_rep*q+BegTrial-1);
%     dft_full_rep{q} = dft_trial_all(trials_per_rep*(q-1)+BegTrial : trials_per_rep*q+BegTrial-1);
%     dft_vel_rep{q} = dft_velocity_all(trials_per_rep*(q-1)+BegTrial : trials_per_rep*q+BegTrial-1);
%     dft_acc_rep{q} = dft_acceleration_all(trials_per_rep*(q-1)+BegTrial : trials_per_rep*q+BegTrial-1);
% end
% 
% % Now construct matrices where each column is a direction and each row is a repetition, then perform the ANOVAs
% % (added the DFT anova matrices so they could be computed within the same loop  -CRF 3/2007)
% resp_mat_anova = [];
% resp_mat_anova_full = [];
% resp_mat_anova_vel = [];
% resp_mat_anova_acc = [];
% for k = 1:length(unique_condition_num)
%     clear select_rep;
%     for q = 1:repetitions
%         n = 0;
%         for i = 1:length(unique_azimuth)
%             for j = 1:length(unique_elevation)          % SELECT now restricted to only 0 deg fixation trials -- CRF 3/2007
%                 select_rep{q} = logical( azimuth_rep{q}==unique_azimuth(i) & elevation_rep{q}==unique_elevation(j) & condition_num_rep{q}==unique_condition_num(k) & fix_x_rep{q}==0 & fix_y_rep{q}==0 );
%                 % length(find(select_rep{q}==1));   % what is this for? -CRF
%                 if sum(select_rep{q}) > 0
%                     n = n+1;
%                     if sum(select_rep{q}) == 1
%                         resp_mat_anova{k}(q,n) = spike_rates_rep{q}(select_rep{q})';
%                         resp_mat_anova_full{k}(q,n) = dft_full_rep{q}(select_rep{q})';   % DFT anova matrices added here
%                         resp_mat_anova_vel{k}(q,n) = dft_vel_rep{q}(select_rep{q})';
%                         resp_mat_anova_acc{k}(q,n) = dft_acc_rep{q}(select_rep{q})';
%                     elseif sum(select_rep{q}) > 1
%                         a = spike_rates_rep{q}(select_rep{q})';%encoutered a cell with 8 repitiions that didnt work
%                         b = dft_full_rep{q}(select_rep{q})'
%                         c = dft_vel_rep{q}(select_rep{q})';
%                         d = dft_acc_rep{q}(select_rep{q})';
%                         resp_mat_anova{k}(q,n) = a(1);
%                         resp_mat_anova_full{k}(q,n) = b(1);
%                         resp_mat_anova_vel{k}(q,n) = c(1);
%                         resp_mat_anova_acc{k}(q,n) = d(1);
%                     end
%                 end
%             end
%         end
%     end
%     [p_anova, table, stats] = anova1(resp_mat_anova{k},[],'off');
%     [p_anova_full, table_full, stats_full] = anova1(resp_mat_anova_full{k},[],'off');
%     [p_anova_vel, table_vel, stats_vel] = anova1(resp_mat_anova_vel{k},[],'off');
%     [p_anova_acc, table_acc, stats_acc] = anova1(resp_mat_anova_acc{k},[],'off');
%     
%     P_anova(k) = p_anova;
%     P_anova_full(k) = p_anova_full;
%     P_anova_vel(k) = p_anova_vel;
%     P_anova_acc(k) = p_anova_acc;
%     
%     anova_table{k} = table;
%     anova_table_full{k} = table_full;
%     anova_table_vel{k} = table_vel;
%     anova_table_acc{k} = table_acc;
%         
%     F_val(k) = anova_table{k}(2,5);
%     F_val_full(k) = anova_table_full{k}(2,5);
%     F_val_vel(k) = anova_table_vel{k}(2,5);
%     F_val_acc(k) = anova_table_acc{k}(2,5);
% end
% 
% F_val = cell2mat(F_val);
% F_val_full = cell2mat(F_val_full);
% F_val_vel = cell2mat(F_val_vel);
% F_val_acc = cell2mat(F_val_acc);


%----------------------------------------------------------------------------
%The Following Anova Code segment, modified by Narayan 02/08/07...
%takes all available trials into account for anova calculation
%----------------------------------------------------------------------------

% Fixed to account for vary_fixation blocks  -CRF 3/2007 
if length(unique_fix_x) >= length(unique_fix_y)
    num_eye_pos = length(unique_fix_x);
else
    num_eye_pos = length(unique_fix_y);
end
trials_per_rep = 26 * length(unique_condition_num) * length(unique_fix_x) * length(unique_fix_y) + num_eye_pos; % the last term corresponds to the null trial(s), one for each eye position
repetitions = ceil( (EndTrial-(BegTrial-1)) / trials_per_rep); %maximum number of repetitions;

%determine the maximum number repetitions over all the trials(directions)..
maxrep = 0;
for k = 1:length(unique_condition_num)
    for i = 1:26 %number of unique directions                % SELECT now restricted to only 0 deg fixation trials -- CRF 3/2007
         select_rep = find(temp_azimuth==unique_azimuth(plot_col(i)) & temp_elevation==unique_elevation(plot_row(i)) & temp_condition_num==unique_condition_num(k) & temp_fix_x==0 & temp_fix_y==0);
         maxrep = max(maxrep, length(select_rep));          % And changed to temp_azimuth, etc. to ensure partial reps are included  -CRF 3/26/2007
    end
end

trial_resp_matrix = [];
trial_full_matrix = [];
trial_vel_matrix = [];
trial_acc_matrix = [];
clear P_anova P_anova_full P_anova_vel P_anova_acc;
clear anova_table anova_table_full anova_table_vel anova_table_acc;
clear F_val F_val_full F_val_vel F_val_acc;

% Fixed to include null trials ('temp' added to vectors in SELECT and spike_rates), for congruency with previous edits  -CRF 3/21/2007
for k=1: length(unique_condition_num)
    for i=1:26 %number of unique directions                             % SELECT now restricted to only 0 deg fixation trials -- CRF 3/2007
         select_rep = find(temp_azimuth==unique_azimuth(plot_col(i)) & temp_elevation==unique_elevation(plot_row(i)) & temp_condition_num==unique_condition_num(k) & temp_fix_x==0 & temp_fix_y==0);
         difc = maxrep - length(select_rep); %difference between actual no. of repetitions and maximum number of repetitions..
         %the above variable is needed to know the number of NaN's that the data needs to be padded with.
         %usually length of select_rep is equal to number of repetitions.. if not pad
         %the rest of the matrix with NaNs for consistent anova analysis
         trial_resp_matrix{k}(:,i) = [temp_spike_rates(select_rep) NaN*ones(1,difc)];
         trial_full_matrix{k}(:,i) = [dft_trial_all(select_rep) NaN*ones(1,difc)];
         trial_vel_matrix{k}(:,i) = [dft_velocity_all(select_rep) NaN*ones(1,difc)];
         trial_acc_matrix{k}(:,i) = [dft_acceleration_all(select_rep) NaN*ones(1,difc)];
         selection{k}(:,i) = [select_rep NaN*ones(1,difc)];
    end
    
    [p_anova, table_resp, stats_resp] = anova1(trial_resp_matrix{k},[],'off');
    [p_anova_full, table_full, stats_full] = anova1(trial_full_matrix{k},[],'off');
    [p_anova_vel, table_vel, stats_vel] = anova1(trial_vel_matrix{k},[],'off');
    [p_anova_acc, table_acc, stats_acc] = anova1(trial_acc_matrix{k},[],'off');

    P_anova(k) = p_anova;                   % Added full DFT anovas, and the F values  -CRF
    P_anova_full(k) = p_anova_full;
    P_anova_vel(k) = p_anova_vel;
    P_anova_acc(k) = p_anova_acc;
    
    anova_table{k} = table_resp;
    anova_table_full{k} = table_full;
    anova_table_vel{k} = table_vel;
    anova_table_acc{k} = table_acc;
        
    F_val(k) = anova_table{k}(2,5);
    F_val_full(k) = anova_table_full{k}(2,5);
    F_val_vel(k) = anova_table_vel{k}(2,5);
    F_val_acc(k) = anova_table_acc{k}(2,5);
    
end

F_val = cell2mat(F_val);
F_val_full = cell2mat(F_val_full);
F_val_vel = cell2mat(F_val_vel);
F_val_acc = cell2mat(F_val_acc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[DDI_mfr, var_term1] = Compute_DDI_anuk(azimuth1, elevation1, spike_rates_dum);  % added by CRF 3/2007
[DDI_full, var_term1]= Compute_DDI_anuk(azimuth1, elevation1, dft_trial);        % added by CRF 3/2007
[DDI_vel, var_term1] = Compute_DDI_anuk(azimuth1, elevation1, dft_velocity);
[DDI_acc, var_term1] = Compute_DDI_anuk(azimuth1, elevation1, dft_acceleration);

% make sure they are the same length before permutation test
if length(dft_velocity) ~= length(spike_rates_dum) | length(dft_velocity) ~= length(dft_trial)
    disp('ERROR - DFT/MFR vectors must be the same= length to compute DDIs and spatial correlation');
    return
end

%permute
for i = 1:1000
    ran_ind = randperm(length(dft_velocity));
    [DDI_mfr_rand(i), var_term1] = Compute_DDI_anuk(azimuth1, elevation1, spike_rates_dum(ran_ind));
    [DDI_full_rand(i), var_term1] = Compute_DDI_anuk(azimuth1, elevation1, dft_acceleration(ran_ind));
    [DDI_vel_rand(i), var_term1] = Compute_DDI_anuk(azimuth1, elevation1, dft_velocity(ran_ind));
    [DDI_acc_rand(i), var_term1] = Compute_DDI_anuk(azimuth1, elevation1, dft_acceleration(ran_ind));
end
t_ddi_mfr = length(find(DDI_mfr_rand>DDI_mfr));
p_ddi_mfr = 2*t_ddi_mfr/1000;
t_ddi_full = length(find(DDI_full_rand>DDI_full));
p_ddi_full = 2*t_ddi_full/1000;
t_ddi_vel = length(find(DDI_vel_rand>DDI_vel));
p_ddi_vel = 2*t_ddi_vel/1000;
t_ddi_acc = length(find(DDI_acc_rand>DDI_acc));
p_ddi_acc = 2*t_ddi_acc/1000;

%% Correlation
resp_mat = [];
resp_mat_vector = [];
resp_mat_std = [];
for i=1:length(unique_azimuth)
    for j=1:length(unique_elevation)
        for k=1
            select = logical( (azimuth1==unique_azimuth(i)) & (elevation1==unique_elevation(j))  );
            if (sum(select) > 0)                
                resp_mat(k, j, i) = mean(dft_velocity(select));
                resp_mat_vector(k, j, i) = mean(dft_velocity(select)); % for vector sum only
                for t = 1 : length(dft_velocity(select));              % this is to calculate response matrix based on each trial
                    spike_temp = dft_velocity(select);                 % in order to calculate error between trials in one condition
                    resp_mat_trial{k}(t, j, i) = spike_temp( t ); % t represents how many repetions each condition
                end
                resp_mat_std(k, j, i) = std(dft_velocity(select));     % calculate std between trials for later DSI usage
            else
                resp_mat(k, j, i) = resp_mat(k,j,1);
                resp_mat_vector(k,j,1) =0; % for vector sum only
                resp_mat_std(k, j, i) = 0;
            end
        end        
    end
end

for i=1:length(unique_azimuth)
    for j=1:length(unique_elevation)
        for k=2
            select = logical( (azimuth1==unique_azimuth(i)) & (elevation1==unique_elevation(j)));
            if (sum(select) > 0)                
                resp_mat(k, j, i) = mean(dft_acceleration(select));
                resp_mat_vector(k, j, i) = mean(dft_acceleration(select)); % for vector sum only
                for t = 1 : length(dft_acceleration(select));              % this is to calculate response matrix based on each trial
                    spike_temp = dft_acceleration(select);                 % in order to calculate error between trials in one condition
                    resp_mat_trial{k}(t, j, i) = spike_temp( t );     % t represents how many repetions each condition
                end
                resp_mat_std(k, j, i) = std(dft_acceleration(select));     % calculate std between trials for later DSI usage
            else                
                resp_mat(k, j, i) = resp_mat(k,j,1);
                resp_mat_vector(k,j,1) =0; % for vector sum only
                resp_mat_std(k, j, i) = 0;
            end
        end        
    end
end

for i=1:length(unique_azimuth)
    for j=1:length(unique_elevation)
        for k=3                   % changed to azimuth1/elevation1 (using spike_rates_dum), for consistency  -CRF
            select = logical( (azimuth1==unique_azimuth(i)) & (elevation1==unique_elevation(j)));
            if (sum(select) > 0)                
                resp_mat(k, j, i) = mean(spike_rates_dum(select));
                resp_mat_vector(k, j, i) = mean(spike_rates_dum(select)); % for vector sum only
                for t = 1 : length(spike_rates_dum(select));              % this is to calculate response matrix based on each trial
                    spike_temp = spike_rates_dum(select);                 % in order to calculate error between trials in one condition
                    resp_mat_trial{k}(t, j, i) = spike_temp( t );     % t represents how many repetions each condition
                end
                resp_mat_std(k, j, i) = std(spike_rates_dum(select));     % calculate std between trials for later DSI usage
            else                
                resp_mat(k, j, i) = resp_mat(k,j,1);
                resp_mat_vector(k,j,1) =0; % for vector sum only
                resp_mat_std(k, j, i) = 0;
            end
        end        
    end
end

%% creat a real 3-D based plot where the center correspond to forward and both lateral edges correspond to backward
%%%% Usually, axis azimuth from left is 270-225-180-135-90--0---90 %%%% 
resp_mat_tran(:,:,1) = resp_mat(:,:,7);
resp_mat_tran(:,:,2) = resp_mat(:,:,6);
resp_mat_tran(:,:,3) = resp_mat(:,:,5);
resp_mat_tran(:,:,4) = resp_mat(:,:,4);
resp_mat_tran(:,:,5) = resp_mat(:,:,3);
resp_mat_tran(:,:,6) = resp_mat(:,:,2);
resp_mat_tran(:,:,7) = resp_mat(:,:,1);
resp_mat_tran(:,:,8) = resp_mat(:,:,8);
resp_mat_tran(:,:,9) = resp_mat_tran(:,:,1);

% flip up-down so that -90 is up and +90 is down (intuitive, and consistent with previous analyses)  -CRF 3/2007
resp_mat_tran_vel = flipud(squeeze(resp_mat_tran(1,:,:)));
resp_mat_tran_acc = flipud(squeeze(resp_mat_tran(2,:,:)));
resp_mat_tran_mfr = flipud(squeeze(resp_mat_tran(3,:,:)));

% % ----------------------------------------------------------------------------
% % TEMPORARY: save figures locally to make sure things make sense  -CRF 3/26/07
%  
% figure; contourf(resp_mat_tran_vel);
% set(gca, 'ytick', [1 2 3 4 5] );
% set(gca, 'ydir' , 'reverse');
% set(gca, 'xticklabel','270|225|180|135|90|45|0|-45|-90');
% set(gca, 'yticklabel','90|45|0|-45|-90');
% title([FILE '  DFTR Vel']);
% colorbar;
% saveas(gcf,['C:\MATLAB6p5\work\Anuk\' FILE '_vel.fig']);
% close;
% 
% figure; contourf(resp_mat_tran_acc);
% set(gca, 'ytick', [1 2 3 4 5] );
% set(gca, 'ydir' , 'reverse');
% set(gca, 'xticklabel','270|225|180|135|90|45|0|-45|-90');
% set(gca, 'yticklabel','90|45|0|-45|-90');
% title([FILE '  DFTR Acc']);
% colorbar;
% saveas(gcf,['C:\MATLAB6p5\work\Anuk\' FILE '_acc.fig']);
% close;
% 
% figure; contourf(resp_mat_tran_mfr);
% set(gca, 'ytick', [1 2 3 4 5] );
% set(gca, 'ydir' , 'reverse');
% set(gca, 'xticklabel','270|225|180|135|90|45|0|-45|-90');
% set(gca, 'yticklabel','90|45|0|-45|-90');
% title([FILE '  MFR']);
% colorbar;
% saveas(gcf,['C:\MATLAB6p5\work\Anuk\' FILE '_mfr.fig']);
% close;
% % ----------------------------------------------------------------------------


%Correlate the mean firing rate contour plot with the dftr contour plot.
count = 1;
for i = 1:length(unique_azimuth)
    for j = 1:length(unique_elevation)
        if (unique_elevation(j) == -90 | unique_elevation(j) == 90) & unique_azimuth(i)~= 0 
            k = 1;
        else
            vel_data(count) = resp_mat_tran(1,j,i);
            acc_data(count) = resp_mat_tran(2,j,i);
            mfr_data(count) = resp_mat_tran(3,j,i);
            count = count+1;
        end
    end
end

vel_mfr = 999;
p_vel_mfr = 999;
acc_mfr = 999;
p_acc_mfr = 999;
acc_vel = 999;
p_acc_vel = 999;

% if P_anova(1)<.05 & P_anova_vel(1)<.05
%     [r,p] = corrcoef( vel_data, mfr_data ); % vel
%     vel_mfr = r(1,2);
%     p_vel_mfr = p(1,2);
% end
% 
% if P_anova(1)<.05 & P_anova_acc(1)<.05
%     [r,p] = corrcoef( acc_data, mfr_data ); % vel
%     acc_mfr = r(1,2);
%     p_acc_mfr = p(1,2);
% end
% 
% if P_anova_vel(1)<.05 & P_anova_acc(1)<.05
%     [r,p] = corrcoef( acc_data, vel_data ); % vel
%     acc_vel = r(1,2);
%     p_acc_vel = p(1,2);
% end

[r,p] = corrcoef( vel_data, mfr_data ); % vel_mfr
vel_mfr = r(1,2);
p_vel_mfr = p(1,2);

[r,p] = corrcoef( acc_data, mfr_data ); % acc_mfr
acc_mfr = r(1,2);
p_acc_mfr = p(1,2);

[r,p] = corrcoef( acc_data, vel_data ); % acc_vel
acc_vel = r(1,2);
p_acc_vel = p(1,2);


% -----------------------------------------------------------------------------------------------------------
% output to text files, 'data1' for each direction, and 'data2' for each cell

% foldername = ('Z:\LabTools\Matlab\TEMPO_Analysis\CommonTools\Frequency_Analysis_data\');
foldername = ('C:\Frequency_Analysis\');

outfile1 = [foldername 'Translation_Frequency_data1.dat'];
outfile2 = [foldername 'Translation_Frequency_data2.dat'];

for i = 1:26
    sprint_txt = ['%s\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t'];
    buff = sprintf(sprint_txt, FILE, unique_azimuth(plot_col(i)), unique_elevation(plot_row(i)), fourier_ratio(i), dft_p(i), fourier_ratio_vel(i), dft_p_vel(i), fourier_ratio_acc(i), dft_p_acc(i), max_phase(i), shift_win(i)*.05);  
    fid = fopen(outfile1, 'a');
    fprintf(fid, '%s', buff);
    fprintf(fid, '\r\n');
    fclose(fid);
end

sprint_txt = ['%s\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t%f\t %f\t %f\t'];
buff = sprintf(sprint_txt, FILE, P_anova(current_stim), P_anova_full(current_stim), P_anova_vel(current_stim), P_anova_acc(current_stim), DDI_mfr, DDI_full, DDI_vel, DDI_acc, p_ddi_mfr, p_ddi_full, p_ddi_vel, p_ddi_acc, vel_mfr, p_vel_mfr,...
    acc_mfr, p_acc_mfr, acc_vel, p_acc_vel, mfr_pref_az, mfr_pref_el, mfr_pref_amp, pref_az_full, pref_el_full, pref_amp_full, pref_az_vel, pref_el_vel, pref_amp_vel, pref_az_acc, pref_el_acc, pref_amp_acc);  
fid = fopen(outfile2, 'a');
fprintf(fid, '%s', buff);
fprintf(fid, '\r\n');
fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Under this, most part is written by Katsu, in order to plot PSTH and
%  Hilbert and DFTR
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% change 26 trajectories to azi x8 and ele x5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

plot_col = [1 1 1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 6 6 6 7 7 7 8 8 8];
plot_row = [5 4 3 2 1 4 3 2 4 3 2 4 3 2 4 3 2 4 3 2 4 3 2 4 3 2];

for i = 1:8
    for j= 1:5
        if (i~=1)
            if (j==1)
                g_dat{i,j}=gauss_dat{5};%j=1,i=2,..8
                Hil{i,j}=Hilbert{5};
            elseif (j==5)
                g_dat{i,j}=gauss_dat{1};%j=5,i=2,..8
                Hil{i,j}=Hilbert{1};
            else
                g_dat{i,j} = gauss_dat{intersect(find(plot_col==i),find(plot_row==j))}; %j=2,3,4,i=2..8
                Hil{i,j}=Hilbert{intersect(find(plot_col==i),find(plot_row==j))};
            end
        else
            g_dat{i,j} = gauss_dat{intersect(find(plot_col==i),find(plot_row==j))}; %i=1, j=1 2 3 4 5
            Hil{i,j}=Hilbert{intersect(find(plot_col==i),find(plot_row==j))};
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pol = -1;%since we want polarity to be positive for phases near 180 and negative for phases near 0
%%% change 26 trajectories to 8 x azi and 5 x ele
for i = 1:8
    for j= 1:5
        if (i~=1)
            if (j==1)
                f_r(i,j)=fourier_ratio(5);%j=1,i=2,..8
                max_p(i,j)=max_phase(5);
%                 vel_f_r(i,j) = pol*fourier_ratio(5)*cos(max_phase(5));
%                 acc_f_r(i,j) = fourier_ratio(5)*sin(max_phase(5));
                vel_f_r(i,j) = fourier_ratio_vel(5);
                acc_f_r(i,j) = fourier_ratio_acc(5);
                d_p(i, j)=dft_p(5);
                d_p_vel(i, j)=dft_p_vel(5);
                d_p_acc(i, j)=dft_p_acc(5);
                s_w(i, j)= shift_win(5);
%                 g_time{i, j}= gausstime{5};
            elseif (j==5)
                f_r(i,j)=fourier_ratio(1);%j=5,i=2,..8
                max_p(i,j)=max_phase(1);
%                 vel_f_r(i,j) = pol*fourier_ratio(1)*cos(max_phase(1));
%                 acc_f_r(i,j) = fourier_ratio(1)*sin(max_phase(1));
                vel_f_r(i,j) = fourier_ratio_vel(1);
                acc_f_r(i,j) = fourier_ratio_acc(1);
                d_p(i, j)=dft_p(1);
                d_p_vel(i, j)=dft_p_vel(1);
                d_p_acc(i, j)=dft_p_acc(1);
                s_w(i, j)= shift_win(1);
%                  g_time{i, j}= gausstime{1};
            else
                f_r(i,j) = fourier_ratio(intersect(find(plot_col==i),find(plot_row==j))); %j=2,3,4,i=2..8
                max_p(i,j) = max_phase(intersect(find(plot_col==i),find(plot_row==j)));
%                 vel_f_r(i,j) = pol*fourier_ratio(intersect(find(plot_col==i),find(plot_row==j)))*cos(max_phase(intersect(find(plot_col==i),find(plot_row==j))));
%                 acc_f_r(i,j) = fourier_ratio(intersect(find(plot_col==i),find(plot_row==j)))*sin(max_phase(intersect(find(plot_col==i),find(plot_row==j))));
                vel_f_r(i,j) = fourier_ratio_vel(intersect(find(plot_col==i),find(plot_row==j)));
                acc_f_r(i,j) = fourier_ratio_acc(intersect(find(plot_col==i),find(plot_row==j)));
                d_p(i, j)=dft_p(intersect(find(plot_col==i),find(plot_row==j)));
                d_p_vel(i, j)=dft_p_vel(intersect(find(plot_col==i),find(plot_row==j)));
                d_p_acc(i, j)=dft_p_acc(intersect(find(plot_col==i),find(plot_row==j)));
                s_w(i, j)= shift_win(intersect(find(plot_col==i),find(plot_row==j)));
%                  g_time{i, j}= gausstime{intersect(find(plot_col==i),find(plot_row==j))};
            end
        else
            f_r(i,j) = fourier_ratio(intersect(find(plot_col==i),find(plot_row==j))); %i=1, j=1 2 3 4 5
            max_p(i,j) = max_phase(intersect(find(plot_col==i),find(plot_row==j)));
%             vel_f_r(i,j) = pol*fourier_ratio(intersect(find(plot_col==i),find(plot_row==j)))*cos(max_phase(intersect(find(plot_col==i),find(plot_row==j))));
%             acc_f_r(i,j) = fourier_ratio(intersect(find(plot_col==i),find(plot_row==j)))*sin(max_phase(intersect(find(plot_col==i),find(plot_row==j))));
            vel_f_r(i,j) = fourier_ratio_vel(intersect(find(plot_col==i),find(plot_row==j)));
            acc_f_r(i,j) = fourier_ratio_acc(intersect(find(plot_col==i),find(plot_row==j)));
            d_p(i, j)=dft_p(intersect(find(plot_col==i),find(plot_row==j)));
            d_p_vel(i, j)=dft_p_vel(intersect(find(plot_col==i),find(plot_row==j)));
            d_p_acc(i, j)=dft_p_acc(intersect(find(plot_col==i),find(plot_row==j)));
            s_w(i, j)= shift_win(intersect(find(plot_col==i),find(plot_row==j)));
%              g_time{i, j}= gausstime{intersect(find(plot_col==i),find(plot_row==j))};
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_start = [StartEventBin(1,1)/timebin, StartEventBin(1,1)/timebin];
x_stop =  [StopEventBin(1,1)/timebin,  StopEventBin(1,1)/timebin];
y_marker=[0,max_count];
% define figure
figure(2);
set(2,'Position', [5,5 1000,700], 'Name', '3D Direction Tuning');
orient landscape;
title(FILE);
axis off;

xoffset=0;
yoffset=0;

% now plot


    axes('position',[0 0 1 1]); 
    xlim([-50,50]);
    ylim([-50,50]);
%     text(-50+xoffset*100,52+yoffset*110, h_title{current_stim} );
    
    text (-50, 44, 'DFTR resp/ p');
    text (-50, 42, 'DFTR vel / p');
    text (-50, 40, 'DFTR acc / p');
    text (-48, 38, 'Delay (sec)');

    text(-49, 35, '-90');text(-49, 20, '-45');text(-49, 5, '0');text(-49, -10, '45');text(-49, -25, '90');
%     
    text(-47,-30, 'Azimth:                 270                     225                    180                      135                    90                       45                       0                        315');
%      
%     values=num2str([P_anova(current_stim), P_anova_vel(current_stim), P_anova_acc(current_stim), mfr_pref_az, mfr_pref_el, mfr_pref_amp, pref_az_vel, pref_el_vel, pref_amp_vel, pref_az_acc, pref_el_acc, pref_amp_acc]);
%     text (-49, -38, 'P-anova-mfr, P-anova-vel, P-anova-acc, mfr-pref-az, mfr-pref-el, mfr-pref-amp, pref-az-vel, pref-el-vel, pref-amp-vel, pref-az-acc, pref-el-acc, pref-amp-acc');
%     text(-49, -42, values)
    text (-40, -33, '            P-anova,     p-DDI,          DDI,          pref-Azimuth,  pref-Elevation,  pref-amplitude');
    text (-48, -36, 'MFR');
    text (-48, -39, 'Fourier Ratio');
    text (-48, -42, 'Velocity(cos)');
    text (-48, -45, 'Accel.(sin)');
    
    values_mfr=num2str([P_anova(current_stim), p_ddi_mfr, DDI_mfr, mfr_pref_az, mfr_pref_el, mfr_pref_amp]);
    values_full=num2str([P_anova_full(current_stim), p_ddi_full, DDI_full, pref_az_full, pref_el_full, pref_amp_full]);
    values_vel=num2str([P_anova_vel(current_stim), p_ddi_vel, DDI_vel, pref_az_vel, pref_el_vel, pref_amp_vel]);
    values_acc=num2str([P_anova_acc(current_stim), p_ddi_acc, DDI_acc, pref_az_acc, pref_el_acc, pref_amp_acc]);
    
    text (-35, -36, values_mfr);
    text (-35, -39, values_full);
    text (-35, -42, values_vel);
    text (-35, -45, values_acc);


    text(25,46, 'Translation');
    axis off;
    hold on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    t_vec =0:0.05:4.95;% by Naryan to plot data x-ordinate
    Hil_vec=1:0.05:2.95;% by Katsu to plot data HilbertTransform
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % add -90 and 90 for j=2,3,4,5,6,7,8
    for i=2:8
        count_y{i, 1, current_stim} = count_y{1, 1, current_stim};
        count_y{i, 5, current_stim} = count_y{1, 5, current_stim};
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % change order 270-225-180-135-90-45-0-315
     for i=1:length(unique_azimuth)                 %
        for j=1:length(unique_elevation)
            if (i < 8 )                                 
               new_count_y{i,j,current_stim}=count_y{8-i,j,current_stim}; 
               new_s_w(i,j)=s_w(8-i,j);
               new_g_dat{i,j}=g_dat{8-i,j};
               new_Hil{i,j}=Hil{8-i,j};
                new_f_r(i,j)=f_r(8-i,j);
                new_max_p(i,j)=max_p(8-i,j);
                new_vel_f_r(i,j)=vel_f_r(8-i,j);
                new_acc_f_r(i,j)=acc_f_r(8-i,j);
                new_d_p(i, j)=d_p(8-i, j);
                new_d_p_vel(i, j)=d_p_vel(8-i, j);
                new_d_p_acc(i, j)=d_p_acc(8-i, j);
            elseif(i==8)
               new_count_y{i,j,current_stim}=count_y{i,j,current_stim}; 
               new_s_w(i,j)=s_w(i,j);
               new_g_dat{i,j}=g_dat{i,j};
               new_Hil{i,j}=Hil{i,j};
                new_f_r(i,j)=f_r(i,j);
                new_max_p(i,j)=max_p(i,j);
                new_vel_f_r(i,j)=vel_f_r(i,j);
                new_acc_f_r(i,j)=acc_f_r(i,j);
                new_d_p(i, j)=d_p(i, j);
                new_d_p_vel(i, j)=d_p_vel(i, j);
                new_d_p_acc(i, j)=d_p_acc(i, j);
            else
               new_count_y{i,j,current_stim}=count_y{7,j,current_stim}; 
               new_s_w(i,j)=s_w(7,j);
               new_g_dat{i,j}=g_dat{7,j};
               new_Hil{i,j}=Hil{7,j};
                new_f_r(i,j)=f_r(7,j);
                new_max_p(i,j)=max_p(7,j);
                new_vel_f_r(i,j)=vel_f_r(7,j);
                new_acc_f_r(i,j)=acc_f_r(7,j);
                new_d_p(i, j)=d_p(7, j);
                new_d_p_vel(i, j)=d_p_vel(7, j);
                new_d_p_acc(i, j)=d_p_acc(7, j);
            end
        end
     end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    for i=1:length(unique_azimuth)                 % 
        for j=1:length(unique_elevation)
%             axes('position',[0.05*i+0.01+xoffset (0.92-0.07*j)+yoffset 0.045 0.045]); 
            axes('position',[0.1*i+0.01+xoffset (0.92-0.14*j)+yoffset 0.09 0.09]); 

bar(t_vec,new_count_y{i,j,current_stim}, 'facecolor', [0.7 0.7 0.7], 'edgecolor', [0.7 0.7 0.7]  );% after rearrangement
            hold on;
%      
% plot( Hil_vec+(new_s_w(i,j)*.05), new_g_dat{i,j}, 'b', 'linewidth',1);
bar( Hil_vec+(new_s_w(i,j)*.05), new_g_dat{i,j}, 'facecolor', [0.5 0.5 0.5], 'edgecolor', [0.5 0.5 0.5]);
            hold on;
plot( Hil_vec, new_Hil{i,j}, 'r', 'linewidth',1.5);
            set( gca, 'xticklabel', ' ' );
            
             plot( x_start*.05, y_marker, 'b-');%This is 0 to 100 should be
%             change to 0 to 5
            plot( x_stop*.05,  y_marker, 'b-');
            
text(2.5, max_count*.8, num2str(new_s_w(i,j)*.05))
            
            % set the same scale for all plot
%             xlim([0,5]);
            xlim([0.5,3.5]);%% recorded 5 sec, middle 2 sec = stimli 1-3 sec (red), +0.5 sec before and after 
            ylim([0,max_count]);
        end    
    end
hold off;
% plot data
    for i=1:length(unique_azimuth)                 % aizmuth 270 are plotted two times in order to make circular data
        for j=1:length(unique_elevation)
%             axes('position',[0.05*i+0.01+xoffset (0.92-0.07*j)+yoffset 0.045 0.045]); 
            axes('position',[0.1*i+0.01+xoffset (0.92-0.14*j)+yoffset+0.09 0.09 0.04]); 
          axis off;
            xlim([0, 100]);
            ylim([0, 100]);
            
            text(0, 100, num2str(new_f_r(i,j)));
            text(0, 60, num2str(new_vel_f_r(i,j)));
            text(0, 20, num2str(new_acc_f_r(i,j)));
            
            text(60, 100, num2str(new_d_p(i, j)));
            text(60, 60, num2str(new_d_p_vel(i, j)));
            text(60, 20, num2str(new_d_p_acc(i, j)));
%             text(10,length(unique_condition_num)-k,'Rotation');
%     text(20,length(unique_condition_num)-k, h_text{k} );
        end    
    end 
%     xoffset=xoffset+0.46;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3);
set(2,'Position', [5,5 1000,700], 'Name', 'DFTR contourf');
orient landscape;
title(FILE);
axis off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear all;
% str1 =  ' %s%f%f%f%f%f%f'; %6 values for direction specific (26/cell for the rotation and direction tuning protocols)
%                            [name  dft_ratio, dft_p, dft_p_vel, dft_p_acc, phase, delay] = textread('Rotation_Frequency_26trajectories.dat', str1);

% str1 =  ' %s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f'; %22 values for cell specific
% [name2 P_anova, P_anova_vel, P_anova_acc, DDI_vel, DDI_acc, p_ddi_vel, p_ddi_acc, corr_vel_mfr, p_vel_mfr, corr_acc_mfr,...
%   p_acc_mfr, corr_acc_vel, p_acc_vel, mfr_pref_az, mfr_pref_el, mfr_pref_amp, pref_az_vel, pref_el_vel, pref_amp_vel,...
%   pref_az_acc, pref_el_acc, pref_amp_acc] = textread('Rotation_Frequency_cell.dat', str1);

%%%%%%%%%%%%%%%%%%5 compare up (following text by Katsu) and down  (This
%%%%%%%%%%%%%%%%%%                                                        file putput part) to convert variables 

%     buff = sprintf(sprint_txt, FILE, fourier_ratio(i), dft_p(i), dft_p_vel(i), dft_p_acc(i), max_phase(i), shift_win(i)*.05); 
%     outfile = [BASE_PATH
%     'ProtocolSpecific\MOOG\rotation3d\Rotation_Frequency_26trajectories.d
%     at'];% one variable includes 26 values

%     buff = sprintf(sprint_txt, FILE, P_anova(current_stim), P_anova_vel(current_stim), P_anova_acc(current_stim), DDI_vel, DDI_acc,...
%     p_ddi_vel, p_ddi_acc, vel_mfr, p_vel_mfr,...
%     acc_mfr, p_acc_mfr, acc_vel, p_acc_vel, mfr_pref_az, mfr_pref_el,
%     mfr_pref_amp, pref_az_vel, pref_el_vel, pref_amp_vel, pref_az_acc, pref_el_acc, pref_amp_acc);  
%      only one value, but P_anova(1) should be to vetibular

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
dft_ratio=fourier_ratio;% 26 values
phase=max_phase;

vel_comp = -dft_ratio.*cos(phase);
acc_comp = dft_ratio.*sin(phase);

%%%%%%%%%%%%%%%% remove {n} !!! only one cell %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   insert all 26 trajectories
% cellno=length(dft_ratio)/26;%n = No. of cell, temporaly 46 cells (03/07/07)
% for n=1:cellno;%n = No. of cell, temporaly 46 cells (03/07/07)
%     for m=1:26;% m=trajectories
%         dft_ratio26(m)=dft_ratio(m);
%         dft_ratio26_vel(m)=vel_comp(m);
%         dft_ratio26_acc(m)=acc_comp(+m);
%         dft_p26(m)=dft_p(m);
%         dft_p26_vel(m)=dft_p_vel(+m);
%         dft_p26_acc(m)=dft_p_acc(m);
%     end
% end % does need for...end?
        dft_ratio26=dft_ratio;
        dft_ratio26_vel=vel_comp;
        dft_ratio26_acc=acc_comp;
        dft_p26=dft_p;
        dft_p26_vel=dft_p_vel;
        dft_p26_acc=dft_p_acc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
count_sig_dft=0;
count_sig_dft_vel=0;
count_sig_dft_acc=0;
% for n=1:(length(dft_ratio)/26);%n = No. of cell, temporaly 46 cells (03/07/07)
    for m=1:26;% m=trajectories
        if dft_p26(m)<0.05
            count_sig_dft=count_sig_dft+1;
        end
        if dft_p26_vel(m)<0.05
            count_sig_dft_vel=count_sig_dft_vel+1;
        end
        if dft_p26_acc(m)<0.05
            count_sig_dft_acc=count_sig_dft_acc+1;
        end
    end
%     count_sig_dft(n)
% end
% P_anova=P_anova(current_stim)
% P_anova ;
% count_sig_dft=count_sig_dft';
% count_sig_dft_vel=count_sig_dft_vel';
% count_sig_dft_acc=count_sig_dft_acc';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% convert i, j
%%% change 26 trajectories to azix8 and elex5
plot_col = [1 1 1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 6 6 6 7 7 7 8 8 8];
plot_row = [5 4 3 2 1 4 3 2 4 3 2 4 3 2 4 3 2 4 3 2 4 3 2 4 3 2];
% for n=1:cellno

for i = 1:8
    for j= 1:5
        if (i~=1)
            if (j==1)
                DFT_mfr(j,i)=dft_ratio26(5);%j=1,i=2,..8 order  j,i =5*8 for contourf
                DFT_vel(j,i)=dft_ratio26_vel(5);
                DFT_acc(j,i)=dft_ratio26_acc(5);
                p_mfr(j,i)=dft_p26(5);
                p_vel(j,i)=dft_p26_vel(5);
                p_acc(j,i)=dft_p26_acc(5);
%                 Hil{i,j}=Hilbert{5};
            elseif (j==5)
                DFT_mfr(j,i)=dft_ratio26(1);%j=5,i=2,..8
                DFT_vel(j,i)=dft_ratio26_vel(1);
                DFT_acc(j,i)=dft_ratio26_acc(1);
                p_mfr(j,i)=dft_p26(1);
                p_vel(j,i)=dft_p26_vel(1);
                p_acc(j,i)=dft_p26_acc(1);
%                 Hil{i,j}=Hilbert{1};
            else
                DFT_mfr(j,i)=dft_ratio26(intersect(find(plot_col==i),find(plot_row==j))); %j=2,3,4,i=2..8
                DFT_vel(j,i)=dft_ratio26_vel(intersect(find(plot_col==i),find(plot_row==j)));
                DFT_acc(j,i)=dft_ratio26_acc(intersect(find(plot_col==i),find(plot_row==j)));
                p_mfr(j,i)=dft_p26(intersect(find(plot_col==i),find(plot_row==j)));
                p_vel(j,i)=dft_p26_vel(intersect(find(plot_col==i),find(plot_row==j)));
                p_acc(j,i)=dft_p26_acc(intersect(find(plot_col==i),find(plot_row==j)));
%                 Hil{i,j}=Hilbert{intersect(find(plot_col==i),find(plot_row==j))};
            end
        else
            DFT_mfr(j,i)=dft_ratio26(intersect(find(plot_col==i),find(plot_row==j))); %i=1, j=1 2 3 4 5
            DFT_vel(j,i)=dft_ratio26_vel(intersect(find(plot_col==i),find(plot_row==j)));
            DFT_acc(j,i)=dft_ratio26_acc(intersect(find(plot_col==i),find(plot_row==j)));
            p_mfr(j,i)=dft_p26(intersect(find(plot_col==i),find(plot_row==j)));
            p_vel(j,i)=dft_p26_vel(intersect(find(plot_col==i),find(plot_row==j)));
            p_acc(j,i)=dft_p26_acc(intersect(find(plot_col==i),find(plot_row==j)));
%             Hil{i,j}=Hilbert{intersect(find(plot_col==i),find(plot_row==j))};
        end
    end
end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for n=1:cellno
for i = 1:8
    for j= 1:5
        if p_mfr(j,i) < 0.05
            p_mfr_p(j,i)=1;
        else
            p_mfr_p(j,i)=0;
        end
        if p_vel(j,i) < 0.05
            p_vel_p(j,i)=1;
        else
            p_vel_p(j,i)=0;
        end
        if p_acc(j,i) < 0.05
            p_acc_p(j,i)=1;
        else
            p_acc_p(j,i)=0;
        end
    end
end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % change order 270-225-180-135-90-45-0-315
     for i=1:length(unique_azimuth)+1                 %
        for j=1:length(unique_elevation)
            if (i < 8 )                                 
               DFT_mfr270(j,i)=DFT_mfr(j,8-i); 
               DFT_vel270(j,i)=DFT_vel(j,8-i);
               DFT_acc270(j,i)=DFT_acc(j,8-i);
               p_mfr_p270(j,i)=p_mfr_p(j,8-i);
               p_vel_p270(j,i)=p_vel_p(j,8-i);
               p_acc_p270(j,i)=p_acc_p(j,8-i);
               
            elseif(i==8)
               DFT_mfr270(j,i)=DFT_mfr(j,i); 
               DFT_vel270(j,i)=DFT_vel(j,i);
               DFT_acc270(j,i)=DFT_acc(j,i);
               p_mfr_p270(j,i)=p_mfr_p(j,i);
               p_vel_p270(j,i)=p_vel_p(j,i);
               p_acc_p270(j,i)=p_acc_p(j,i);
%                 new_d_p_acc(i, j)=d_p_acc(i, j);
            else
               DFT_mfr270(j,i)=DFT_mfr(j,7); 
               DFT_vel270(j,i)=DFT_vel(j,7);
               DFT_acc270(j,i)=DFT_acc(j,7);
               p_mfr_p270(j,i)=p_mfr_p(j,7);
               p_vel_p270(j,i)=p_vel_p(j,7);
               p_acc_p270(j,i)=p_acc_p(j,7);
            end
        end
     end
%%%%%   plot  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for f=1:ceil(cellno./4);% ceil (46 cells/4) is 11.5 turned to 12 figures
%     figure(f+5);orient landscape;
%     for p=1:4;
        
            subplot(3,2,1);contourf( DFT_mfr270(:,:) );colorbar;set(gca, 'ydir' , 'reverse');
                 t=['Spatial Tuning of response / p-full = ',num2str( P_anova_full(current_stim) )];title(t)
                  set(gca, 'XTickMode','manual');
                    set(gca, 'xtick',[1,2,3,4,5,6,7,8,9]);
                    set(gca, 'xticklabel','270|225|180|135|90|45|0|-45|-90'); 
                    xlabel('Azimuth');
                  set(gca, 'YTickMode','manual');
                    set(gca, 'ytick',[1,2,3,4,5]);
                    set(gca, 'yticklabel','-90|-45|0|45|90');
                    ylabel('Elevation');
            subplot(3,2,2);contourf( p_mfr_p270(:,:) );colorbar('YTickLabel',{'NotSig',' ',' ',' ',' ',' '});set(gca, 'ydir' , 'reverse');
                t=['# of Sig Directions = ',num2str( count_sig_dft )];title(t);
                  set(gca, 'XTickMode','manual');
                    set(gca, 'xtick',[1,2,3,4,5,6,7,8,9]);
                    set(gca, 'xticklabel','270|225|180|135|90|45|0|-45|-90'); 
                    xlabel('Azimuth');
                  set(gca, 'YTickMode','manual');
                    set(gca, 'ytick',[1,2,3,4,5]);
                    set(gca, 'yticklabel','-90|-45|0|45|90');
                    ylabel('Elevation');
%                 set(gca, 'ytick', [] ); 
            subplot(3,2,3);contourf( DFT_vel270(:,:) );colorbar;set(gca, 'ydir' , 'reverse');
                t=['Spatial Tuning of Velocity / p-Vel = ',num2str( P_anova_vel(current_stim))];title(t)
                set(gca, 'XTickMode','manual');
                    set(gca, 'xtick',[1,2,3,4,5,6,7,8,9]);
                    set(gca, 'xticklabel','270|225|180|135|90|45|0|-45|-90'); 
                    xlabel('Azimuth');
                  set(gca, 'YTickMode','manual');
                    set(gca, 'ytick',[1,2,3,4,5]);
                    set(gca, 'yticklabel','-90|-45|0|45|90');
                    ylabel('Elevation');
%                 set(gca, 'ytick', [] ); 
            subplot(3,2,4);contourf( p_vel_p270(:,:) );colorbar('YTickLabel',{'NotSig',' ',' ',' ',' ',' '});set(gca, 'ydir' , 'reverse');
                t=['# of Sig Directions = ',num2str( count_sig_dft_vel )];title(t);
                set(gca, 'XTickMode','manual');
                    set(gca, 'xtick',[1,2,3,4,5,6,7,8,9]);
                    set(gca, 'xticklabel','270|225|180|135|90|45|0|-45|-90'); 
                    xlabel('Azimuth');
                  set(gca, 'YTickMode','manual');
                    set(gca, 'ytick',[1,2,3,4,5]);
                    set(gca, 'yticklabel','-90|-45|0|45|90');
                    ylabel('Elevation');
%                 set(gca, 'ytick', [] ); 
            subplot(3,2,5);contourf( DFT_acc270(:,:) );colorbar;set(gca, 'ydir' , 'reverse');
                t=['Spatial Tuning of Acceleration / p-Acc = ',num2str( P_anova_acc(current_stim) )];title(t)
                set(gca, 'XTickMode','manual');
                    set(gca, 'xtick',[1,2,3,4,5,6,7,8,9]);
                    set(gca, 'xticklabel','270|225|180|135|90|45|0|-45|-90'); 
                    xlabel('Azimuth');
                  set(gca, 'YTickMode','manual');
                    set(gca, 'ytick',[1,2,3,4,5]);
                    set(gca, 'yticklabel','-90|-45|0|45|90');
                    ylabel('Elevation');
%                 set(gca, 'ytick', [] ); 
            subplot(3,2,6);contourf( p_acc_p270(:,:) );colorbar('YTickLabel',{'NotSig',' ',' ',' ',' ',' '});set(gca, 'ydir' , 'reverse');
                t=['# of Sig Directions = ',num2str( count_sig_dft_acc )];title(t);
                set(gca, 'XTickMode','manual');
                    set(gca, 'xtick',[1,2,3,4,5,6,7,8,9]);
                    set(gca, 'xticklabel','270|225|180|135|90|45|0|-45|-90'); 
                    xlabel('Azimuth');
                  set(gca, 'YTickMode','manual');
                    set(gca, 'ytick',[1,2,3,4,5]);
                    set(gca, 'yticklabel','-90|-45|0|45|90');
                    ylabel('Elevation');
%                 set(gca, 'ytick', [] ); 
   axes('position',[0 0 1 1]); 
    xlim([-50,50]);
    ylim([-50,50]);
    text(5,47, 'Translation');
    text(-10,47,FILE);
    axis off;
    hold on;
%     end
% end

% 
% pause % to print out during batch files
return;
