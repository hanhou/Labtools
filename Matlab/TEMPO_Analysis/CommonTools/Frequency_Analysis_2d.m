%----------------------------------------------------------------------------------------------------------------------
%Frequency_Analysis.m calculates the DFTR, phase, latency, velocity and acceleration
%components, preferred directions (for MFR, velocity and acceleration),
%DDI's, vector sum amplitudes, correlation of 3d tuning and anova for neuronal responses IN THE
%VESITIBULAR (stim_type = 1) CONDITION.

function Frequency_Analysis_2d(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
Path_Defs;
ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP

%get the column of values for azimuth and elevation and stim_type
temp_azimuth = data.moog_params(AZIMUTH,:,MOOG);
temp_elevation = data.moog_params(ELEVATION,:,MOOG);
temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG); 
temp_amplitude = data.moog_params(AMPLITUDE,:,MOOG); 
temp_spike_data = data.spike_data(1,:);
dummy_spike_data = temp_spike_data;
%mean firing rate of each trial depending on the start and stop offsets
temp_spike_rates = data.spike_rates(SpikeChan, :);

%get indices of any NULL conditions (for measuring spontaneous activity
null_trials = logical( (temp_azimuth == data.one_time_params(NULL_VALUE)) );

%now, remove trials from direction and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:length(temp_azimuth);		%a vector of trial indices
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
spike_rates= temp_spike_rates(~null_trials & select_trials);

%baseline firing rate
spon_found = find(null_trials==1);
spon_resp = mean(temp_spike_rates(null_trials));
for i = 1:length(spon_found)
    total_null(i) = temp_spike_rates(spon_found(i));
end
spon_count = floor(sum(temp_spike_rates(null_trials))/20);

unique_azimuth = munique(azimuth');
unique_elevation = munique(elevation');
unique_stim_type = munique(stim_type');

condition_num = stim_type;
h_title{1}='Vestibular';
h_title{2}='Visual';
h_title{3}='Combined';
unique_condition_num = munique(condition_num');

timebin=50; % Bins the mfr for that 50msec. So 40 bins total for the 2 seconds of folded data.
% sample frequency depends on test duration
frequency=length(temp_spike_data)/length(select_trials); % actually its 5000 now for 5 sec 
% length of x-axis
x_length = frequency/timebin;
% x-axis for plot PSTH
x_time=1:(frequency/timebin);

total_trials = floor( length(azimuth) /(26*length(unique_condition_num)) );
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
time_step=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
current_stim = 1;%vestibular=1, visual= 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k=1: length(unique_condition_num)  % it's the stim_type
        for i=1: length(unique_azimuth)
            select = logical( (azimuth==unique_azimuth(i)) & (elevation==0) & (condition_num==unique_condition_num(k)) );            
            % get rid off -90 and 90 cases
            if (sum(select) > 0)
                resp{k}(i) = mean(spike_rates(select)); 
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
                        count_condition{i,repeat} = dummy_count{repeat};
                    end
                end
                
                count_y_trial{i,k}(:,:) = temp_count;  % each trial's PSTH in vestibular condition
                if k == 1
                    count_y_trial1{i}(:,:) = temp_count;
                end
                
                dim=size(temp_count);
                if dim(1) > 1;
                    count_y{i,k} = mean(temp_count);
                else
                    count_y{i,k}= temp_count;     % for only one repetition cases
                end
            else
                resp{k}(i) = 0; 
                count_y{i,k}=0;
            end   
            % normalize count_y
            if max(count_y{i,k})~=0;
                count_y_norm{i,k}=count_y{i,k} / max(count_y{i,k});
            else
                count_y_norm{i,k}=0;
            end
        end
    % now find the peak
    row_max = find( resp{k}==max(resp{k}) );
    % it is likely there are two peaks with same magnitude, choose the first one arbitraly
    row_m{k}=row_max(1);
    if max(count_y{row_max(1), k})~=0;
        count_y_max{k} = count_y{row_max(1), k} / max(count_y{row_max(1), k});
    else
        count_y_max{k} =0;
    end
    % find the largest y to set scale later
    if max(count_y{row_max(1), k}) > max_count
        max_count = max(count_y{row_max(1), k});
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%row_m{1};
%col_m{1};

%find empty trials
sum_dat = zeros(1,100);
for j = 1:length(unique_azimuth)
    sum_dat = sum_dat + count_y{j,current_stim};
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
else%this is one of katsu's cells
    span = 80;
end

%perform frequency analysis
for j=1:length(unique_azimuth)
    gdat = count_y{j,current_stim};
    temp_dat{j} = gdat;
    temp_dat1{j} = gdat;
    %---------------------------------------
    %find hilbert transform to get peak.
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
    dummy_dat{j} = hildata(18:63) - mean(hildata(14:20));
    total_data{j} = gdat(1:span);        
    hil{j} = hilbert(dummy_dat{j});% find the peak within the window that we are considering (1.5 sec)
    abb = abs(hil{j});
    hhh{j} = abb;
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
    
    % stim peak time = 2 sec = 40th bin
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
for i = 1:length(unique_azimuth)
    DFT_ratio(i) = fourier_ratio(i);
end
% %--------------------------------------------------------------------------
%calculate preferred direction for velocity and acceleration by resolving
%the phase into two components - 180/0 and -90/90
vel_sum = zeros(1,8);
acc_sum = zeros(1,8);
for i = 1:length(unique_azimuth)
    pol = -1;%since we want polarity to be positive for phases near 180 and negative for phases near 0
    temp_phase = max_phase(i);          
    vel_sum(i) = pol*fourier_ratio(i)*cos(temp_phase);
    acc_sum(i) = fourier_ratio(i)*sin(temp_phase);
end

[pref_az_vel pref_el_vel pref_amp_vel] = vectorsum(vel_sum);
[pref_az_acc pref_el_acc pref_amp_acc] = vectorsum(acc_sum);
[mfr_pref_az mfr_pref_el mfr_pref_amp] = vectorsum(resp{current_stim});
%------------------------------------------------

%find the dft ratio of each individual trial
trial_time = 0:.05:5;
Azims = [0 45 90 135 180 225 270 315];
Elevs = [-90 -45 0 45 90];
azimuth1 = azimuth(find(condition_num == 1));
elevation1 = elevation(find(condition_num == 1));
spike_rates_dum = spike_rates(find(condition_num == 1));
temp_spike_data = data.spike_data(1,:);

ii = 1;
for i = 1:5000:length(temp_spike_data)
   tdat = temp_spike_data(i:i+4999);
   step = 1:50:5001;
   for n=1:(x_length)
       all_count{ii}(n) = sum(tdat(step(n):step(n+1)-1));
   end
   ii = ii+1;
end

for k=1:length(unique_condition_num)
        for i=1: length(unique_azimuth)
            select = logical( (temp_azimuth==unique_azimuth(i)) & (temp_elevation==0) & (temp_stim_type==k) );%Changed it to temp_azimuth as the variable azimuth doesnt have null trials in it.
            %replacements azimuth<=>temp_azimuth, elevation <=>
            %temp_elevation, condition_num <=> temp_stim_type
            % get rid off -90 and 90 cases
            if (sum(select) > 0)
                act_found = find( select==1 );
                for repeat=1:length(act_found) 
                    temp_trial_data = all_count{act_found(repeat)};%changes for each trial
                    %trial_peak = peak_time(find( (Azims(plot_col) == unique_azimuth(i)) & (Elevs(plot_row) == 0) )); % same for each trial
                    trial_span = span;% same for each trial
                    signal = temp_trial_data(14:20);%changes for each trial
                    rand_sig = signal(ran_in{i});% same for each trial
                    shift_trial = shift_win(i);% same for each trial
                    temp_trial_data = temp_trial_data(1:trial_span);
                    
                    if shift_trial > 8 & span == 68 % if we shift by more than 8 bins, there is not recorded data at the end of the window
                        for ii = 1:(shift_trial - 8)
                            temp_trial_data(trial_span+ii) = rand_sig(ii);
                        end
                    end
                    %                 now we have to calculate the dft ratio for each trial
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

for i=1: length(unique_azimuth)
        select = logical( (azimuth1==unique_azimuth(i)) & (elevation1==0) );            
        % get rid off -90 and 90 cases
        if (sum(select) > 0)
            act_found = find( select==1 );
            for repeat=1:length(act_found) 
                temp_trial_data = count_y_trial1{i}(repeat,:);%changes for each trial
                %trial_peak = peak_time(find( (Azims(plot_col) == unique_azimuth(i)) & (Elevs(plot_row) == unique_elevation(j)) )); % same for each trial
                trial_span = span;% same for each trial
                signal = temp_trial_data(14:20);%changes for each trial
                rand_sig = signal(ran_in{i});% same for each trial
                shift_trial = shift_win(i);% same for each trial
                temp_trial_data = temp_trial_data(1:trial_span);
                
                if shift_trial > 8 & span == 68 % if we shift by more than 8 bins, there is not recorded data at the end of the window
                    for ii = 1:(shift_trial - 8)
                        temp_trial_data(trial_span+ii) = rand_sig(ii);
                    end
                end
                %                 now we have to calculate the dft ratio for each trial
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

% % % 
%  %ANOVA modified by Aihua, it does not require whole trials, it does not matter if trial stopped during repetition
trials_per_rep = 8 * length(unique_condition_num) + 1;%1 denotes the null trial.
repetitions = floor( (EndTrial-(BegTrial-1)) / trials_per_rep);
% first parse raw data into repetitions, including null trials
for q = 1:repetitions
    azimuth_rep{q} = temp_azimuth(trials_per_rep*(q-1)+BegTrial : trials_per_rep*q+BegTrial-1);
    elevation_rep{q} = temp_elevation(trials_per_rep*(q-1)+BegTrial : trials_per_rep*q+BegTrial-1);
    condition_num_rep{q} = temp_stim_type(trials_per_rep*(q-1)+BegTrial : trials_per_rep*q+BegTrial-1);
    spike_rates_rep{q} = temp_spike_rates(trials_per_rep*(q-1)+BegTrial : trials_per_rep*q+BegTrial-1);
    dft_vel_rep{q} = dft_velocity_all(trials_per_rep*(q-1)+BegTrial : trials_per_rep*q+BegTrial-1);
    dft_acc_rep{q} = dft_acceleration_all(trials_per_rep*(q-1)+BegTrial : trials_per_rep*q+BegTrial-1);
end
resp_mat_anova = [];
for k=1: length(unique_condition_num)
    clear select_rep;
    for q=1:1:repetitions
        n = 0;
        for i=1:length(unique_azimuth)
            select_rep{q} = logical( azimuth_rep{q}==unique_azimuth(i) & elevation_rep{q}==0 & condition_num_rep{q}==unique_condition_num(k) );
            length(find(select_rep{q}==1));
            if (sum(select_rep{q}) > 0)
                n = n+1;
                if sum(select_rep{q}) == 1
                    resp_mat_anova{k}(q,n) = spike_rates_rep{q}(select_rep{q})';
                elseif sum(select_rep{q}) > 1
                    a = spike_rates_rep{q}(select_rep{q})';%encoutered a cell with 8 repitiions that didnt work
                    resp_mat_anova{k}(q,n) = a(1);
                end
            end
        end
    end
    [p_anova, table, stats] = anova1(resp_mat_anova{k},[],'off');
    P_anova(k) = p_anova;
    anova_table{k} = table;
    F_val(k) = anova_table{k}(2,5);
end
F_val = cell2mat(F_val);
% % 
% % 
resp_mat_anova_vel = [];
resp_mat_anova_acc = [];
for k=1: length(unique_condition_num)
    clear select_rep;
    for q=1:1:repetitions
        n = 0;
        for i=1:length(unique_azimuth)
                select_rep{q} = logical( azimuth_rep{q}==unique_azimuth(i) & elevation_rep{q}==0 & condition_num_rep{q}==unique_condition_num(k) );
                length(find(select_rep{q}==1));
                if (sum(select_rep{q}) > 0)
                    n = n+1;
                    if sum(select_rep{q}) == 1
                        resp_mat_anova_vel{k}(q,n) = dft_vel_rep{q}(select_rep{q})';
                        resp_mat_anova_acc{k}(q,n) = dft_acc_rep{q}(select_rep{q})';
                    elseif sum(select_rep{q}) > 1
                        a = dft_vel_rep{q}(select_rep{q})';%encoutered a cell with 8 repitiions that didnt work
                        b = dft_acc_rep{q}(select_rep{q})';
                        resp_mat_anova_vel{k}(q,n) = a(1);
                        resp_mat_anova_acc{k}(q,n) = b(1);
                    end
                end
         end
    end
    [p_anova_vel, table_vel, stats_vel] = anova1(resp_mat_anova_vel{k},[],'off');
    [p_anova_acc, table_acc, stats_acc] = anova1(resp_mat_anova_acc{k},[],'off');
    P_anova_vel(k) = p_anova_vel;
    anova_table_vel{k} = table_vel;
    F_val_vel(k) = anova_table_vel{k}(2,5);
    P_anova_acc(k) = p_anova_acc;
    anova_table_acc{k} = table_acc;
    F_val_acc(k) = anova_table_acc{k}(2,5);
end
F_val_vel = cell2mat(F_val_vel);
F_val_acc = cell2mat(F_val_acc);

[DDI_vel, var_term1] = Compute_DDI_2d(azimuth1, elevation1, dft_velocity);
[DDI_acc, var_term1] = Compute_DDI_2d(azimuth1, elevation1, dft_acceleration);

%permute
for i = 1:1000
    ran_ind = randperm(length(dft_velocity));
    [DDI_vel_rand(i), var_term1] = Compute_DDI_2d(azimuth1, elevation1, dft_velocity(ran_ind));
    [DDI_acc_rand(i), var_term1] = Compute_DDI_2d(azimuth1, elevation1, dft_acceleration(ran_ind));
end
t_ddi_vel = length(find(DDI_vel_rand>DDI_vel));
p_ddi_vel = 2*t_ddi_vel/1000;
t_ddi_acc = length(find(DDI_acc_rand>DDI_acc));
p_ddi_acc = 2*t_ddi_acc/1000;

%% Correlation
resp_mat = [];
for i=1:length(unique_azimuth)
        for k=1
            select = logical( (azimuth1==unique_azimuth(i)) & (elevation1==0)  );
            if (sum(select) > 0)                
                resp_mat(k, j, i) = mean(dft_velocity(select));
                resp_mat_vector(k, j, i) = mean(dft_velocity(select)); % for vector sum only
                for t = 1 : length(dft_velocity(select));              % this is to calculate response matrix based on each trial
                    spike_temp = dft_velocity(select);                 % in order to calculate error between trials in one condition
                    resp_mat_trial{k}(t, j, i) = spike_temp( t );     % t represents how many repetions each condition
                end
                resp_mat_std(k, j, i) = std(dft_velocity(select));     % calculate std between trials for later DSI usage
            else
                resp_mat(k, j, i) = resp_mat(k,j,1);
                resp_mat_vector(k,j,1) =0; % for vector sum only
                resp_mat_std(k, j, i) = 0;
            end
        end        
end

for i=1:length(unique_azimuth)
        for k=2
            select = logical( (azimuth1==unique_azimuth(i)) & (elevation1==0));
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

for i=1:length(unique_azimuth)
        for k=3
            select = logical( (azimuth==unique_azimuth(i)) & (elevation==0) & (condition_num == 1 ));
            if (sum(select) > 0)                
                resp_mat(k, j, i) = mean(spike_rates(select));
                resp_mat_vector(k, j, i) = mean(spike_rates(select)); % for vector sum only
                for t = 1 : length(spike_rates(select));              % this is to calculate response matrix based on each trial
                    spike_temp = spike_rates(select);                 % in order to calculate error between trials in one condition
                    resp_mat_trial{k}(t, j, i) = spike_temp( t );     % t represents how many repetions each condition
                end
                resp_mat_std(k, j, i) = std(spike_rates(select));     % calculate std between trials for later DSI usage
            else
                resp_mat(k, j, i) = resp_mat(k,j,1);
                resp_mat_vector(k,j,1) =0; % for vector sum only
                resp_mat_std(k, j, i) = 0;
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

%Correlate the mean firing rate contour plot with the dftr contour plot.

for i = 1:length(unique_azimuth)
            vel_data(i) = resp_mat_tran(1,j,i);
            acc_data(i) = resp_mat_tran(2,j,i);
            mfr_data(i) = resp_mat_tran(3,j,i);
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
[r,p] = corrcoef( vel_data, mfr_data ); % vel
vel_mfr = r(1,2);
p_vel_mfr = p(1,2);

[r,p] = corrcoef( acc_data, mfr_data ); % vel
acc_mfr = r(1,2);
p_acc_mfr = p(1,2);

[r,p] = corrcoef( acc_data, vel_data ); % vel
acc_vel = r(1,2);
p_acc_vel = p(1,2);

foldername = ('Z:\Users\Ayanna\frequency_analysis\');

for i = 1:8
    sprint_txt = ['%s\t        %f\t        %f\t       %f\t      %f\t        %f\t        %f\t'];
    buff = sprintf(sprint_txt, FILE, fourier_ratio(i), dft_p(i), dft_p_vel(i), dft_p_acc(i), max_phase(i), shift_win(j)*.05);  
    outfile = [foldername '1DAzimuth_Frequency_data1.dat'];
    fid = fopen(outfile, 'a');
    fprintf(fid, '%s', buff);
    fprintf(fid, '\r\n');
    fclose(fid);
end

sprint_txt = ['%s\t        %f\t        %f\t       %f\t      %f\t        %f\t        %f\t        %f\t        %f\t       %f\t      %f\t        %f\t        %f\t       %f\t        %f\t       %f\t      %f\t        %f\t        %f\t        %f\t        %f\t       %f\t      %f\t'];
buff = sprintf(sprint_txt, FILE, P_anova(1), P_anova_vel(1), P_anova_acc(1), DDI_vel, DDI_acc, p_ddi_vel, p_ddi_acc, vel_mfr, p_vel_mfr,...
    acc_mfr, p_acc_mfr, acc_vel, p_acc_vel, mfr_pref_az, mfr_pref_el, mfr_pref_amp, pref_az_vel, pref_el_vel, pref_amp_vel, pref_az_acc, pref_el_acc, pref_amp_acc);  
outfile = [foldername '1DAzimuth_Frequency_data2.dat'];
fid = fopen(outfile, 'a');
fprintf(fid, '%s', buff);
fprintf(fid, '\r\n');
fclose(fid);
return;


% str1 =  ' %s%f%f%f%f%f%f';
% [name  dft_ratio, dft_p, dft_p_vel, dft_p_acc, phase, delay] = textread('Frequency_data1.dat', str1);
% str1 =  ' %s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f';
% [name2 P_anova, P_anova_vel, P_anova_acc, DDI_vel, DDI_acc, p_ddi_vel, p_ddi_acc, vel_mfr, p_vel_mfr, corr_acc_mfr, p_acc_mfr, corr_acc_vel, p_acc_vel, corr_mfr_pref_az, mfr_pref_el, mfr_pref_amp, pref_az_vel, pref_el_vel, pref_amp_vel, pref_az_acc, pref_el_acc, pref_amp_acc] = textread('Frequency_data2.dat', str1);
%