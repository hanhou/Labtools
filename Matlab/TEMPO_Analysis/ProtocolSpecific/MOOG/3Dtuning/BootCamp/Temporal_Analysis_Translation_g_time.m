%----------------------------------------------------------------------------------------------------------------------
%Frequency_Analysis.m calculates the DFTR, phase, latency, velocity and acceleration
%components, preferred directions (for MFR, velocity and acceleration),
%DDI's, vector sum amplitudes, correlation of 3d tuning and anova for neuronal responses IN THE
%VESITIBULAR (stim_type = 1) CONDITION.

function Temporal_Analysis_Rotation(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
Path_Defs;
ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP

%get the column of values for azimuth and elevation and stim_type
temp_azimuth = data.moog_params(AZIMUTH,:,MOOG);
temp_elevation = data.moog_params(ELEVATION,:,MOOG);
temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG); 
temp_amplitude = data.moog_params(AMPLITUDE,:,MOOG); 
temp_spike_data = data.spike_data(SpikeChan,:);

%mean firing rate of each trial depending on the start and stop offsets
%temp_spike_rates = data.spike_rates(SpikeChan, :);

% temp_azimuth = data.moog_params(ROT_AZIMUTH,:,MOOG);
% temp_elevation = data.moog_params(ROT_ELEVATION,:,MOOG);
% temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG); 
% temp_spike_data = data.spike_data(SpikeChan,:);

temp_spike_rates = data.spike_rates(SpikeChan, :);    
dummy_spike_data = temp_spike_data;

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

timebin=50;
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

for k=1: length(current_stim)  % it's the stim_type
    for j=1:length(unique_elevation)
        for i=1: length(unique_azimuth)
            select = logical( (azimuth==unique_azimuth(i)) & (elevation==unique_elevation(j)) & (condition_num==unique_condition_num(k)) );            
            % get rid off -90 and 90 cases
            if (sum(select) > 0)
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
row_m{1};
col_m{1};

% order in which the directions are plotted
plot_col = [1 1 1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 6 6 6 7 7 7 8 8 8];
plot_row = [5 4 3 2 1 4 3 2 4 3 2 4 3 2 4 3 2 4 3 2 4 3 2 4 3 2];

%find empty trials
sum_dat = zeros(1,100);
for j = 1:26
    sum_dat = sum_dat + count_y{plot_col(j), plot_row(j),current_stim};
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
for j=1:26
    gdat = count_y{plot_col(j), plot_row(j),current_stim};
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     gaustime=[x_time(20+shift_win(j)):x_time(59+shift_win(j))]
%     gausstime=0.05*gaustime;    
%     %calculate DFT ratio
%     %40 pt DFT with mean removed
%     [f, amp, resp_phase] = FT(gausstime, gauss_dat{j}, 40, 1, 0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5    
    gaustime{j}=[x_time(20+shift_win(j)):x_time(59+shift_win(j))]
    gausstime{j}=0.05*gaustime{j};    
    %calculate DFT ratio
    %40 pt DFT with mean removed
    [f, amp, resp_phase] = FT(gausstime{j}, gauss_dat{j}, 40, 1, 0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
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
% %--------------------------------------------------------------------------
%calculate preferred direction for velocity and acceleration by resolving
%the phase into two components - 180/0 and -90/90
vel_sum = zeros(5,8);
acc_sum = zeros(5,8);
for i = 1:26
    pol = -1;%since we want polarity to be positive for phases near 180 and negative for phases near 0
    temp_phase = max_phase(i);          
    vel_sum(plot_row(i), plot_col(i)) = pol*fourier_ratio(i)*cos(temp_phase);
    acc_sum(plot_row(i), plot_col(i)) = fourier_ratio(i)*sin(temp_phase);
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
    for j=1:length(unique_elevation)
        for i=1: length(unique_azimuth)
            select = logical( (temp_azimuth==unique_azimuth(i)) & (temp_elevation==unique_elevation(j)) & (temp_stim_type==k) );
            %we use temp_azimuth, temp_elevation, temp_stim_type instead of
            %azimuth, elevations and condition_num which do not have any
            %null trials in them.
            % get rid off -90 and 90 cases
            if (sum(select) > 0)
                act_found = find( select==1 );
                for repeat=1:length(act_found) 
                    temp_trial_data = all_count{act_found(repeat)};%changes for each trial
                    trial_peak = peak_time(find( (Azims(plot_col) == unique_azimuth(i)) & (Elevs(plot_row) == unique_elevation(j)) )); % same for each trial
                    trial_span = span;% same for each trial
                    signal = temp_trial_data(14:20);%changes for each trial
                    rand_sig = signal(ran_in{find( (Azims(plot_col) == unique_azimuth(i)) & (Elevs(plot_row) == unique_elevation(j)) )});% same for each trial
                    shift_trial = shift_win(find( (Azims(plot_col) == unique_azimuth(i)) & (Elevs(plot_row) == unique_elevation(j)) ));% same for each trial
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
end

for j=1:length(unique_elevation)
    for i=1: length(unique_azimuth)
        select = logical( (azimuth1==unique_azimuth(i)) & (elevation1==unique_elevation(j)) );            
        % get rid off -90 and 90 cases
        if (sum(select) > 0)
            act_found = find( select==1 );
            for repeat=1:length(act_found) 
                temp_trial_data = count_y_trial1{i,j}(repeat,:);%changes for each trial
                trial_peak = peak_time(find( (Azims(plot_col) == unique_azimuth(i)) & (Elevs(plot_row) == unique_elevation(j)) )); % same for each trial
                trial_span = span;% same for each trial
                signal = temp_trial_data(14:20);%changes for each trial
                rand_sig = signal(ran_in{find( (Azims(plot_col) == unique_azimuth(i)) & (Elevs(plot_row) == unique_elevation(j)) )});% same for each trial
                shift_trial = shift_win(find( (Azims(plot_col) == unique_azimuth(i)) & (Elevs(plot_row) == unique_elevation(j)) ));% same for each trial
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
end



%The Following Anova Code segment, modified by Narayan 02/08/07...
%takes all available trials into account for anova calculation
trials_per_rep = 26 * length(unique_condition_num) + 1;% 1 corresponds to the null trial
repetitions = ceil( (EndTrial-(BegTrial-1)) / trials_per_rep); %maximum number of repetitions;

plot_col = [1 1 1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 6 6 6 7 7 7 8 8 8];
plot_row = [5 4 3 2 1 4 3 2 4 3 2 4 3 2 4 3 2 4 3 2 4 3 2 4 3 2];

%determine the maximum number repetitions over all the trials(directions)..
maxrep = 0;
for k=1: length(unique_condition_num)
    for i=1:26%number of unique directions
         select_rep = find(azimuth==unique_azimuth(plot_col(i)) & elevation==unique_elevation(plot_row(i)) & condition_num==unique_condition_num(k));
         maxrep = max(maxrep, length(select_rep));
     end
 end

for k=1: length(unique_condition_num)
    for i=1:26%number of unique directions
         select_rep = find(azimuth==unique_azimuth(plot_col(i)) & elevation==unique_elevation(plot_row(i)) & condition_num==unique_condition_num(k));
         difc = maxrep - length(select_rep); %difference between actual no. of repetitions and maximum number of repetitions..
         %the above variable is needed to know the number of NaN's that the data needs to be padded with.
         %usually length of select_rep is equal to number of repetitions.. if not pad
         %the rest of the matrix with NaNs for consistent anova analysis
         trial_resp_matrix{k}(:,i) = [spike_rates(select_rep) NaN*ones(1,difc)];
         trial_vel_matrix{k}(:,i)=[dft_velocity_all(select_rep) NaN*ones(1,difc)];
         trial_acc_matrix{k}(:,i)=[dft_acceleration_all(select_rep) NaN*ones(1,difc)];
         selection{k}(:,i) = [select_rep NaN*ones(1,difc)];
    end
    
 [p_anova, table_resp, stats_resp] = anova1(trial_resp_matrix{k},[],'off');
 [p_anova_vel, table_vel, stats_vel] = anova1(trial_vel_matrix{k},[],'off');
 [p_anova_acc, table_acc, stats_acc] = anova1(trial_acc_matrix{k},[],'off');
 P_anova(k) = p_anova;
 P_anova_vel(k) = p_anova_vel;
 P_anova_acc(k) = p_anova_acc;
end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[DDI_vel, var_term1] = Compute_DDI_anuk(azimuth1, elevation1, dft_velocity);
[DDI_acc, var_term1] = Compute_DDI_anuk(azimuth1, elevation1, dft_acceleration);

%permute
for i = 1:1000
    ran_ind = randperm(length(dft_velocity));
    [DDI_vel_rand(i), var_term1] = Compute_DDI_anuk(azimuth1, elevation1, dft_velocity(ran_ind));
    [DDI_acc_rand(i), var_term1] = Compute_DDI_anuk(azimuth1, elevation1, dft_acceleration(ran_ind));
end
t_ddi_vel = length(find(DDI_vel_rand>DDI_vel));
p_ddi_vel = 2*t_ddi_vel/1000;
t_ddi_acc = length(find(DDI_acc_rand>DDI_acc));
p_ddi_acc = 2*t_ddi_acc/1000;

%% Correlation
resp_mat = [];
for i=1:length(unique_azimuth)
    for j=1:length(unique_elevation)
        for k=1
            select = logical( (azimuth1==unique_azimuth(i)) & (elevation1==unique_elevation(j))  );
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
        for k=3
            select = logical( (azimuth==unique_azimuth(i)) & (elevation==unique_elevation(j)) & (condition_num == 1 ));
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
count = 1;
for i = 1:length(unique_azimuth)
    for j = 1:length(unique_elevation)
        if (unique_elevation(j) == -90 | unique_elevation(j) == 90) & unique_azimuth(i)~= 0 
            k =1;
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
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%26 back to azi and ele 8*5=40 trajectories for PSTH plot
% g_dat{1,5}=gauss_dat(1);
% g_dat{1,4}=gauss_dat(2);
% g_dat{1,3}=gauss_dat(3);
% g_dat{1,2}=gauss_dat(4);
% g_dat{1,1}=gauss_dat(5);
% for z=2:8
%     g_dat{z,4}=gauss_dat(z*3);
%     g_dat{z,3}=gauss_dat(z*3+1);
%     g_dat{z,2}=gauss_dat(z*3+2);
%     g_dat{z,5}=gauss_dat(1);
%     g_dat{z,1}=gauss_dat(5);
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% change 26 to azi and ele
plot_col = [1 1 1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 6 6 6 7 7 7 8 8 8];
plot_row = [5 4 3 2 1 4 3 2 4 3 2 4 3 2 4 3 2 4 3 2 4 3 2 4 3 2];

for i = 1:8
    for j= 1:5
        if (i~=1)
            if (j==1)
                g_dat{i,j}=gauss_dat{5};%j=1,i=2,..8
            elseif (j==5)
                g_dat{i,j}=gauss_dat{1};%j=5,i=2,..8
            else
                g_dat{i,j} = gauss_dat{intersect(find(plot_col==i),find(plot_row==j))}; %j=2,3,4,i=2..8
            end
        else
            g_dat{i,j} = gauss_dat{intersect(find(plot_col==i),find(plot_row==j))}; %i=1, j=1 2 3 4 5
        end
    end
end

for i = 1:8
    for j= 1:5
        if (i~=1)
            if (j==1)
                f_r(i,j)=fourier_ratio(5);%j=1,i=2,..8
                d_p(i, j)=dft_p(5);
                d_p_vel(i, j)=dft_p_vel(5);
                d_p_acc(i, j)=dft_p_acc(5);
                s_w(i, j)= shift_win(5);
                g_time{i, j}= gausstime{5};
            elseif (j==5)
                f_r(i,j)=fourier_ratio(1);%j=5,i=2,..8
                d_p(i, j)=dft_p(1);
                d_p_vel(i, j)=dft_p_vel(1);
                d_p_acc(i, j)=dft_p_acc(1);
                s_w(i, j)= shift_win(1);
                 g_time{i, j}= gausstime{1};
            else
                f_r(i,j) = fourier_ratio(intersect(find(plot_col==i),find(plot_row==j))); %j=2,3,4,i=2..8
                d_p(i, j)=dft_p(intersect(find(plot_col==i),find(plot_row==j)));
                d_p_vel(i, j)=dft_p_vel(intersect(find(plot_col==i),find(plot_row==j)));
                d_p_acc(i, j)=dft_p_acc(intersect(find(plot_col==i),find(plot_row==j)));
                s_w(i, j)= shift_win(intersect(find(plot_col==i),find(plot_row==j)));
                 g_time{i, j}= gausstime{intersect(find(plot_col==i),find(plot_row==j))};
            end
        else
            f_r(i,j) = fourier_ratio(intersect(find(plot_col==i),find(plot_row==j))); %i=1, j=1 2 3 4 5
            d_p(i, j)=dft_p(intersect(find(plot_col==i),find(plot_row==j)));
            d_p_vel(i, j)=dft_p_vel(intersect(find(plot_col==i),find(plot_row==j)));
            d_p_acc(i, j)=dft_p_acc(intersect(find(plot_col==i),find(plot_row==j)));
            s_w(i, j)= shift_win(intersect(find(plot_col==i),find(plot_row==j)));
             g_time{i, j}= gausstime{intersect(find(plot_col==i),find(plot_row==j))};
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

    
%     if( xoffset > 0.5)          % now temperarily 2 pictures one row and 2 one column
%         yoffset = yoffset-0.42;
%         xoffset = 0;
%     end
    % output some text 
    axes('position',[0 0 1 1]); 
    xlim([-50,50]);
    ylim([-50,50]);
%     text(-50+xoffset*100,52+yoffset*110, h_title{current_stim} );
    text (-47, 46, 'DFTR');
    text (-47, 44, 'dft p');
    text (-47, 42, 'dft p vel');
    text (-47, 40, 'dft p acc');

%     text(-47,-40, 'Azim: 270       225       180        135        90        45        0        315        270');
    text(-47,-32, 'Azimth:                 270                     225                  180                    135                   90                     45                     0                      315');
%      text(-47,-40, 'Azimth: 0       45       90        135        180       225        270        315        ');  
    values=num2str([P_anova(current_stim), P_anova_vel(current_stim), P_anova_acc(current_stim), mfr_pref_az, mfr_pref_el, mfr_pref_amp, pref_az_vel, pref_el_vel, pref_amp_vel, pref_az_acc, pref_el_acc, pref_amp_acc]);
    text (-49, -38, 'P-anova-mfr, P-anova-vel, P-anova-acc, mfr-pref-az, mfr-pref-el, mfr-pref-amp, pref-az-vel, pref-el-vel, pref-amp-vel, pref-az-acc, pref-el-acc, pref-amp-acc');
    text(-49, -42, values)
    text(25, 46, 'Translation');
    axis off;
    hold on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    t_vec =0:0.05:4.95;% by Naryan to plot data x-ordinate
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % add -90 and 90 for j=2,3,4,5,6,7,8
    for i=2:8
        count_y{i, 1, current_stim} = count_y{1, 1, current_stim};
        count_y{i, 5, current_stim} = count_y{1, 5, current_stim};
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for i=1:length(unique_azimuth)                 % aizmuth 270 are plotted two times in order to make circular data NO
        for j=1:length(unique_elevation)
%             axes('position',[0.05*i+0.01+xoffset (0.92-0.07*j)+yoffset 0.045 0.045]); 
            axes('position',[0.1*i+0.01+xoffset (0.92-0.14*j)+yoffset 0.09 0.09]); 
%            bar(t_vec,count_y{i,j,current_stim});% start with 0 45 90 135
            if (i < 8 )                                 % temporarilly line output figure with contour one, so that the middle panel corresponds to 90 deg,                             
                bar( t_vec,count_y{8-i,j,current_stim} );    % which is forward motion and the lateral edges correspond to 270 deg which is backward motion
            elseif(i==8)
                bar( t_vec,count_y{i,j,current_stim}); 
            else
                bar( t_vec,count_y{7,j,current_stim} ); 
            end
            hold on;
%             plot( x_start, y_marker, 'r-');%This is 0 to 100 should be
%             change to 0 to 5
%             plot( x_stop,  y_marker, 'r-');
            
%             plot( gausstime+(s_w(i,j)*.05), g_dat{i,j}, 'b');%start with 0 45 90 135
%             if (i<8)
%                 plot( gausstime+(s_w(8-i,j)*.05), g_dat{8-i,j}, 'b');%from 270-225-180
%             elseif (i==8)
%                 plot( gausstime+(s_w(i,j)*.05), g_dat{i,j}, 'b');
%             else
%                 plot( gausstime+(s_w(7,j)*.05), g_dat{7,j}, 'b');
%             end

%                if (i<8)
%                 plot( gausstime, g_dat{8-i,j}, 'b');%from 270-225-180    
%             elseif (i==8)
%                 plot( gausstime, g_dat{i,j}, 'b');
%             else
%                 plot( gausstime, g_dat{7,j}, 'b');
%             end

            if (i<8)
                plot( g_time{i, j}, g_dat{8-i,j}, 'b');%from 270-225-180
            elseif (i==8)
                plot( g_time{i, j}, g_dat{i,j}, 'b');
            else
                plot( g_time{i, j}, g_dat{7,j}, 'b');
            end
            hold on;
            set( gca, 'xticklabel', ' ' );
            
             plot( x_start*.05, y_marker, 'r-');%This is 0 to 100 should be
%             change to 0 to 5
            plot( x_stop*.05,  y_marker, 'r-');
            
            % set the same scale for all plot
%             xlim([0,5]);
            xlim([0.5,3.5]);% recorded 5 sec, middle 2 sec = stimli 1-3 sec (red), +0.5 sec before and after 
            ylim([0,max_count]);
        end    
    end
hold off;
    for i=1:length(unique_azimuth)                 % aizmuth 270 are plotted two times in order to make circular data
        for j=1:length(unique_elevation)
%             axes('position',[0.05*i+0.01+xoffset (0.92-0.07*j)+yoffset 0.045 0.045]); 
            axes('position',[0.1*i+0.01+xoffset (0.92-0.14*j)+yoffset+0.09 0.09 0.04]); 
          axis off;
            xlim([0, 100]);
            ylim([0, 100]);
            
            text(30, 100, num2str(f_r(i,j)));
            text(30, 60, num2str(d_p(i, j)));
            text(30, 30, num2str(d_p_vel(i, j)));
            text(30, 00, num2str(d_p_acc(i, j)));
            
            text(70, 45, num2str(s_w(i,j)*.05));
%             text(10,length(unique_condition_num)-k,'Rotation');
%     text(20,length(unique_condition_num)-k, h_text{k} );
        end    
    end 
%     xoffset=xoffset+0.46;




% % foldername = 'Z:\LabTools\Matlab\TEMPO_Analysis\ProtocolSpecific\MOOG\3Dtuning';
% 
% for i = 1:26
%     sprint_txt = ['%s\t        %f\t        %f\t       %f\t      %f\t        %f\t        %f\t'];
%     buff = sprintf(sprint_txt, FILE, fourier_ratio(i), dft_p(i), dft_p_vel(i), dft_p_acc(i), max_phase(i), shift_win(j)*.05);  
%     outfile = [BASE_PATH 'ProtocolSpecific\MOOG\3dtuning\Translation_Frequency_26trajectories.dat'];%Ve
% %     outfile = [BASE_PATH 'ProtocolSpecific\MOOG\3dtuning\Visu_Translation_Frequency_26trajectories.dat'];%Visu
%     fid = fopen(outfile, 'a');
%     fprintf(fid, '%s', buff);
%     fprintf(fid, '\r\n');
%     fclose(fid);
% end
% 
% sprint_txt = ['%s\t        %f\t        %f\t       %f\t      %f\t        %f\t        %f\t        %f\t        %f\t       %f\t      %f\t        %f\t        %f\t       %f\t        %f\t       %f\t      %f\t        %f\t        %f\t        %f\t        %f\t       %f\t      %f\t'];
% buff = sprintf(sprint_txt, FILE, P_anova(current_stim), P_anova_vel(current_stim), P_anova_acc(current_stim), DDI_vel, DDI_acc, p_ddi_vel, p_ddi_acc, vel_mfr, p_vel_mfr,...
%     acc_mfr, p_acc_mfr, acc_vel, p_acc_vel, mfr_pref_az, mfr_pref_el, mfr_pref_amp, pref_az_vel, pref_el_vel, pref_amp_vel, pref_az_acc, pref_el_acc, pref_amp_acc);  
% outfile = [BASE_PATH 'ProtocolSpecific\MOOG\3dtuning\Translation_Frequency_cell.dat'];%vestib
% % outfile = [BASE_PATH 'ProtocolSpecific\MOOG\3dtuning\Visu_Translation_Frequency_cell.dat'];%Visu
% fid = fopen(outfile, 'a');
% fprintf(fid, '%s', buff);
% fprintf(fid, '\r\n');
% fclose(fid);


return;


% str1 =  ' %s%f%f%f%f%f%f';
% [name  dft_ratio, dft_p, dft_p_vel, dft_p_acc, phase, delay] = textread('Frequency_data1.dat', str1);
% str1 =  ' %s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f';
% [name2 P_anova, P_anova_vel, P_anova_acc, DDI_vel, DDI_acc, p_ddi_vel, p_ddi_acc, vel_mfr, p_vel_mfr, corr_acc_mfr, p_acc_mfr, corr_acc_vel, p_acc_vel, corr_mfr_pref_az, mfr_pref_el, mfr_pref_amp, pref_az_vel, pref_el_vel, pref_amp_vel, pref_az_acc, pref_el_acc, pref_amp_acc] = textread('Frequency_data2.dat', str1);
%