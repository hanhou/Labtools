function Lambert_DirectionTuningPlot_3D_fix(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

Path_Defs;
ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP

%get the column of values for azimuth and elevation and stim_type
temp_azimuth = data.moog_params(AZIMUTH,:,MOOG);
temp_elevation = data.moog_params(ELEVATION,:,MOOG);
temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG); 
temp_amplitude = data.moog_params(AMPLITUDE,:,MOOG); 
temp_fix_x    =  data.moog_params(FIX_X,:,MOOG);
temp_fix_y    =  data.moog_params(FIX_Y,:,MOOG);
temp_spike_data = data.spike_data(1,:);
temp_spike_rates = data.spike_rates(SpikeChan, :);
dummy_spike_data = temp_spike_data;

%get indices of any NULL conditions (for measuring spontaneous activity
null_trials = logical( (temp_azimuth == data.one_time_params(NULL_VALUE)) );

%now, remove trials from direction and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:length(temp_azimuth);		% a vector of trial indices
bad_trials = find(temp_spike_rates > 3000);   % cut off 3k frequency which definately is not cell's firing response
%at any given point, cell response cannot be greater than 1
abnormal = find(temp_spike_data > 1);
temp_spike_data(1,abnormal) = 1;

if ( bad_trials ~= NaN)
    select_trials= ( (trials >= BegTrial) & (trials <= EndTrial) & (trials~=bad_trials) );
else 
    select_trials= ( (trials >= BegTrial) & (trials <= EndTrial) ); 
end

azimuth = temp_azimuth(~null_trials & select_trials);
elevation = temp_elevation(~null_trials & select_trials);
stim_type = temp_stim_type(~null_trials & select_trials);
amplitude = temp_amplitude(~null_trials & select_trials);
fix_x     = temp_fix_x(~null_trials & select_trials);
fix_y     = temp_fix_y(~null_trials & select_trials);
spike_rates= temp_spike_rates(~null_trials & (trials >= BegTrial) & (trials <= EndTrial));
% notice that this bad_trials is the number without spon trials 
bad_trials = find(spike_rates > 3000);   % cut off 3k frequency which definately is not cell's firing response

%baseline firing rate 
spon_found = find(null_trials == 1);
spon_resp = mean(temp_spike_rates(null_trials));
spon_count = floor(sum(temp_spike_rates(null_trials))/20);
for i = 1:length(spon_found)
    total_null(i) = temp_spike_rates(spon_found(i));
end

unique_azimuth  = munique(azimuth');
unique_elevation = munique(elevation');
unique_stim_type = munique(stim_type');
unique_amplitude = munique(amplitude');
unique_fix_x    =  munique(fix_x');
unique_fix_y    =  munique(fix_y');

total_trials = floor( length(azimuth) /(length(unique_stim_type)*78) );

if length(unique_fix_y) == 1
    condition_num = fix_x;
    temp_condition_num = temp_fix_x;
else
    condition_num = fix_y; 
    temp_condition_num = temp_fix_y;
end

unique_condition_num = munique(condition_num');

% azimuth = azimuth(find(condition_num==0));
% elevation = elevation(find(condition_num==0));
% stim_type = stim_type(find(condition_num==0));
% amplitude = amplitude(find(condition_num==0));
% spike_rates = spike_rates(find(condition_num==0));

% add parameters here
% timebin for plot PSTH
timebin=50;
% sample frequency depends on test duration
frequency=length(temp_spike_data)/length(select_trials);  
% length of x-axis
x_length = frequency/timebin;
% x-axis for plot PSTH
x_time=1:(frequency/timebin);

stim_duration = length(temp_spike_data)/length(temp_azimuth);                  

% take spontaneous activity out of the whole spike_data
Discard_trials = find(null_trials==1 | trials <BegTrial | trials >EndTrial);
for i = 1 : length(Discard_trials)
    temp_spike_data( 1, ((Discard_trials(i)-1)*stim_duration+1) :  Discard_trials(i)*stim_duration ) = 9999;
end
spike_data = temp_spike_data( temp_spike_data~=9999 );
% count spikes from raster data (spike_data)
max_count = 1;
time_step=1;

% if unique_stim_type(1) == 2
%     stim = 1;
% end
% if length(unique_stim_type) > 1
%     stim = 2;
% end

stim = 1;
if unique_stim_type(1) == 1
    
    for k=1:length(unique_condition_num)   
        for j=1:length(unique_elevation)
            for i=1: length(unique_azimuth)
                select = logical( (azimuth==unique_azimuth(i)) & (elevation==unique_elevation(j)) & (stim_type == 1)  & condition_num ==unique_condition_num(k));            
                % get rid off -90 and 90 cases
                if (sum(select) > 0)
                    resp{k}(j,i) = mean(spike_rates(select));
                    act_found = find( select==1 );
                    % count spikes per timebin on every same condition trials
                    for repeat=1:length(act_found) 
                        for n=1:(x_length)
                            temp_count(repeat,n)=sum(spike_data(1,(frequency*(act_found(repeat)-1)+time_step):(frequency*(act_found(repeat)-1)+n*timebin)));
                            dummy_count{repeat}(n) = sum(spike_data(1,(frequency*(act_found(repeat)-1)+time_step):(frequency*(act_found(repeat)-1)+n*timebin)));
                            time_step=time_step+timebin;
                            
                        end
                        time_step=1;     
                        
                    end
                    
                    if k == 2
                        count_y_trial{i,j}(:,:) = temp_count;  % each trial's PSTH 
                        % get the average of the total same conditions if repetion is > 1
                        % if (length(act_found) > 1);
                    end               
                    
                    dim=size(temp_count);
                    if dim(1) > 1;
                        count_y{i,j,k} = mean(temp_count);
                        %the total number of spikes in the 5 trials
                        %                     total_spike_count{i,j,k} = temp_count;
                    else
                        count_y{i,j,k}= temp_count;     % for only one repetition cases
                        %                     total_spike_count{i,j,k} = temp_count;
                    end
                    %                 if k == 2
                    %                     total_visual_spikes{i,j} = temp_count;
                    %                 end
                    
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
    
    
    
    
    
    
    stim = 1;
    
    % order in which the directions are plotted
    plot_col = [1 1 1 1 1 2 2 2 2 2 3 3 3 3 3 4 4 4 4 4 5 5 5 5 5 6 6 6 6 6 7 7 7 7 7 8 8 8 8 8];
    plot_row = [5 4 3 2 1 5 4 3 2 1 5 4 3 2 1 5 4 3 2 1 5 4 3 2 1 5 4 3 2 1 5 4 3 2 1 5 4 3 2 1];
    
    sum_dat = zeros(1,100);
    is_empty = zeros(1,40);%Add by Cah 04-20-06
    for j = 1:40
        if sum(count_y{plot_col(j), plot_row(j),stim}) == 0
            is_empty(j) = 1;
        end
        if is_empty(j) == 1
            count_y{plot_col(j), plot_row(j),stim} = zeros(1, 100);
        end
        sum_dat = sum_dat + count_y{plot_col(j), plot_row(j),stim};
    end
    
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
    
    for j=1:40
        gdat = count_y{plot_col(j), plot_row(j),stim};
        temp_dat{j} = gdat;
        
        %---------------------------------------
        %find hilbert transform to get peak.
        shift_win(j) = 0;
        peak_time(j) = 40;
        flagged_dir(j) = 0;
        if is_empty(j) == 0
            
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
            dummy_dat{j} = hildata(18:63) - mean(hildata(14:20));%only consider 1.8 to 2.7 sec for finding the peak of the response
            total_data{j} = gdat(1:span);        
            hil = hilbert(dummy_dat{j});% find the peak within the window that we are considering (1.5 sec)
            abb = abs(hil);
            HT = abs(hil) - mean(abb(4:5));% remove the mean to analyze just the peak or else the center of mass will be biased
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
            
            % stim peak time = 2 sec = 40
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
        else
            total_data{j} = gdat;%it doesnt really matter how long this 
        end
        
        for k=1:40
            gaussdata(k)=total_data{j}(k+20+shift_win(j)); %window is shifted depending on how delayed the peak of the envelope (from the hilbert trans) is
        end
        gauss_dat{j} = gaussdata;
        gaustime=[x_time(20+shift_win(j)):x_time(60+shift_win(j))];
        gausstime=0.05*gaustime;    
        %calculate DFT ratio
        if is_empty(j) == 0
            [f, amp, resp_phase] = FT(gausstime, gauss_dat{j}, 40, 1, 0);
            f1 = mean(amp(2:4));
            f2 = mean(amp(5:end));
            if f2 == 0
                f_rat = 0;
            else
                f_rat = f1/f2;
            end
            max_phase(j) = mean(resp_phase(find(amp == max(amp(2:4)))));
            dft_p(j) = DFT_perm(gauss_dat{j}, f_rat, 1000);
        else
            f = 0;
            amp = 0;
            f_rat = 0;
            max_phase(j) = 999;
            dft_p(j) = 999;
        end
        fourier_ratio(j) = f_rat;
        mean_rates(j) = mean(temp_dat{j}(30:50));
    end
    
    
    goo = find(dft_p~=999);
    Min_resp = min(fourier_ratio(goo));
    Max_resp = max(fourier_ratio(goo));
    
    temp_max = find(fourier_ratio == max(fourier_ratio));
    shift_phase = max_phase; % shift negative angle by 360
    angle_add = zeros(1,length(shift_phase));
    angle_add(find(shift_phase<0)) = 2*pi;
    shift_phase = shift_phase+angle_add;
    m_phase = shift_phase(temp_max);
    % calculate the preferred direction using vector sum
    DFT_ratio = zeros(5,8);
    DFT_dummy = zeros(5,8);
    delta_phase = zeros(5,8);
    for i = 1:40
        if sum(fourier_ratio(i)) > 0
            DFT_ratio(plot_row(i), plot_col(i)) = fourier_ratio(i);
            DFT_dummy(plot_row(i), plot_col(i)) = fourier_ratio(i);
            cell_response{plot_row(i), plot_col(i)} = dummy_dat{i};
            delta_phase(plot_row(i), plot_col(i)) = min([abs((m_phase - shift_phase(i))) (2*pi - abs((m_phase - shift_phase(i))))])*180/pi;
        end
        if plot_row(i) == 1 & plot_col(i) ~= 1
            delta_phase(plot_row(i),plot_col(i)) = delta_phase(1,1);
            DFT_dummy(plot_row(i),plot_col(i)) = DFT_ratio(1,1);
        end
        if plot_row(i) == 5 & plot_col(i) ~= 1
            delta_phase(plot_row(i),plot_col(i)) = delta_phase(5,1);
            DFT_dummy(plot_row(i),plot_col(i)) = DFT_ratio(5,1);
        end
    end
    %--------------------------------------------------------------------------
    %calculate preferred direction for velocity and acceleration by resolving
    %the phase into two components - 180/0 and -90/90
    vel_sum = zeros(5,8);
    acc_sum = zeros(5,8);
    for i = 1:40
        pol = -1;%since we want polarity to be positive for phases near 180 and negative for phases near 0
        if is_empty(i) == 0
            temp_phase = max_phase(i);          
            vel_sum(plot_row(i), plot_col(i)) = pol*fourier_ratio(i)*cos(temp_phase);
            acc_sum(plot_row(i), plot_col(i)) = fourier_ratio(i)*sin(temp_phase);
        end
    end
    [pref_az_vel pref_el_vel pref_amp_vel] = vectorsum(vel_sum);
    [pref_az_acc pref_el_acc pref_amp_acc] = vectorsum(acc_sum);
    [pref_az_mfr pref_el_mfr mfr_pref_amp] = vectorsum(resp{stim});
    
    stimu = 1;
    
    %find the DFT ratio for the null trials
    null_index = find(null_trials == 0);
    for i = 1 : length(null_index)
        dummy_spike_data( 1, ((null_index(i)-1)*stim_duration+1) :  null_index(i)*stim_duration ) = 9999;
    end
    null_data = dummy_spike_data( dummy_spike_data~=9999 );
    num_null = length(null_data)/5000;
    %bin the null trials now
    k=0;
    bin_size = 1:50:5000;
    for i = 1:num_null
        null_dat{i} = null_data(k+1: k+ 5000);
        k = k+5000;
        for j = 1:length(bin_size)-1
            binned_null{i}(j) = sum(null_dat{i}(bin_size(j):bin_size(j+1)-1));
        end
        [null_f null_amp nul_phase] = FT(0:.05:1.95, binned_null{i}(21:60), 40, 1, 0);
        f1_null = mean(null_amp(2:4));
        f2_null = mean(null_amp(5:end));
        %     avoid dividing by 0
        if f2_null == 0
            null_ratio(i) = 0;
        else
            null_ratio(i) = f1_null/f2_null;
        end
        null_phase(i) = mean(nul_phase(find(null_amp == max(null_amp(2:4)))));
        null_velocity(i) = -null_ratio(i)*cos(null_phase(i));
        null_acceleration(i) = null_ratio(i)*sin(null_phase(i));
    end
    spon_ratio_vel = mean(null_velocity);
    spon_ratio_acc = mean(null_acceleration);
    
    %find the dft ratio of each individual trial
    trial_time = 0:.05:5;
    Azims = [0 45 90 135 180 225 270 315];
    Elevs = [-90 -45 0 45 90];
    azimuth1 = azimuth(find(stim_type == 1 & condition_num == 0));
    elevation1 = elevation(find(stim_type == 1 & condition_num == 0));
    spike_rates_dum = spike_rates(find(stim_type == 1 & condition_num == 0));
    for j=1:length(unique_elevation)
        for i=1: length(unique_azimuth)
            select = logical( (azimuth1==unique_azimuth(i)) & (elevation1==unique_elevation(j)) );            
            % get rid off -90 and 90 cases
            if (sum(select) > 0)
                act_found = find( select==1 );
                for repeat=1:length(act_found) 
                    temp_trial_data = count_y_trial{i,j}(repeat,:);%changes for each trial
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
                    %                     dft_velocity_const_phase(act_found(repeat)) = -dft_trial(act_found(repeat))*cos(max_phase(find( (Azims(plot_col) == unique_azimuth(i)) & (Elevs(plot_row) == unique_elevation(j))) ));%use amplitude but not phase
                    %                     dft_acceleration_const_phase(act_found(repeat)) = dft_trial(act_found(repeat))*sin(max_phase(find( (Azims(plot_col) == unique_azimuth(i)) & (Elevs(plot_row) == unique_elevation(j))) ));%use amplitude but not phase
                    
                end
            end
        end
    end
    
    
    
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
    dft_velocity_all = zeros(1,length(temp_azimuth));
    dft_acceleration_all = zeros(1,length(temp_azimuth));
    trial_phase_all = zeros(1,length(temp_azimuth));
    
    for k=1:length(unique_stim_type)
        for j=1:length(unique_elevation)
            for i=1: length(unique_azimuth)
                select = logical( (temp_azimuth==unique_azimuth(i)) & (temp_elevation==unique_elevation(j)) & (temp_stim_type==unique_stim_type(k)) );            
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
                        test(act_found(repeat)) = 1;
                    end
                end
            end
        end
    end
    
    
    %  %ANOVA modified by Aihua, it does not require whole trials, it does not matter if trial stopped during repetition
    trials_per_rep = 26 * length(unique_stim_type) + 1;
    repetitions = floor( (EndTrial-(BegTrial-1)) / trials_per_rep);
    % first parse raw data into repetitions, including null trials
    for q = 1:repetitions
        azimuth_rep{q} = temp_azimuth(trials_per_rep*(q-1)+BegTrial : trials_per_rep*q+BegTrial-1);
        elevation_rep{q} = temp_elevation(trials_per_rep*(q-1)+BegTrial : trials_per_rep*q+BegTrial-1);
        stim_type_rep{q} = temp_stim_type(trials_per_rep*(q-1)+BegTrial : trials_per_rep*q+BegTrial-1);
        spike_rates_rep{q} = temp_spike_rates(trials_per_rep*(q-1)+BegTrial : trials_per_rep*q+BegTrial-1);
        dft_vel_rep{q} = dft_velocity_all(trials_per_rep*(q-1)+BegTrial : trials_per_rep*q+BegTrial-1);
        dft_acc_rep{q} = dft_acceleration_all(trials_per_rep*(q-1)+BegTrial : trials_per_rep*q+BegTrial-1);
        condition_num_rep{q} = temp_condition_num(trials_per_rep*(q-1)+BegTrial : trials_per_rep*q+BegTrial-1);
    end
    
    resp_mat_anova = [];
    for k=1: length(unique_stim_type)
        clear select_rep;
        for q=1:1:repetitions
            n = 0;
            for i=1:length(unique_azimuth)
                for j=1:length(unique_elevation)
                    select_rep{q} = logical( azimuth_rep{q}==unique_azimuth(i) & elevation_rep{q}==unique_elevation(j) & stim_type_rep{q}==unique_stim_type(k) & condition_num_rep{q} == 0);
                    if (sum(select_rep{q}) > 0)
                        n = n+1;
                        a = spike_rates_rep{q}(select_rep{q})';
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
    
    
    
    resp_mat_anova_vel = [];
    resp_mat_anova_acc = [];
    for k=1: length(unique_stim_type)
        clear select_rep;
        for q=1:1:repetitions
            n = 0;
            for i=1:length(unique_azimuth)
                for j=1:length(unique_elevation)
                    select_rep{q} = logical( azimuth_rep{q}==unique_azimuth(i) & elevation_rep{q}==unique_elevation(j) & stim_type_rep{q}==unique_stim_type(k) & condition_num_rep{q}==0);
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
    
    %% ADD CODE HERE FOR PLOTTING
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
                    %                 resp_mat_ste(k, j, i) = resp_mat_std(k, j, i)/ sqrt(length(find( (azimuth1==unique_azimuth(i)) & (elevation1==unique_elevation(j))&)) );
                else
                    %                resp_mat_trial{k}(t, j, i) = 0;
                    resp_mat(k, j, i) = resp_mat(k,j,1);
                    resp_mat_vector(k,j,1) =0; % for vector sum only
                    resp_mat_std(k, j, i) = 0;
                    %                 resp_mat_ste(k, j, i) = 0;
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
                    %                 resp_mat_ste(k, j, i) = resp_mat_std(k, j, i)/ sqrt(length(find( (azimuth==unique_azimuth(i)) & (elevation==unique_elevation(j))&(condition_num==1) )) );
                else
                    %                resp_mat_trial{k}(t, j, i) = 0;
                    resp_mat(k, j, i) = resp_mat(k,j,1);
                    resp_mat_vector(k,j,1) =0; % for vector sum only
                    resp_mat_std(k, j, i) = 0;
                    %                 resp_mat_ste(k, j, i) = 0;
                end
            end        
        end
    end
    
    for i=1:length(unique_azimuth)
        for j=1:length(unique_elevation)
            for k=3
                select = logical( (azimuth==unique_azimuth(i)) & (elevation==unique_elevation(j)) & (stim_type==1 ) & (condition_num == 0));
                if (sum(select) > 0)                
                    resp_mat(k, j, i) = mean(spike_rates(select));
                    resp_mat_vector(k, j, i) = mean(spike_rates(select)); % for vector sum only
                    for t = 1 : length(spike_rates(select));              % this is to calculate response matrix based on each trial
                        spike_temp = spike_rates(select);                 % in order to calculate error between trials in one condition
                        resp_mat_trial{k}(t, j, i) = spike_temp( t );     % t represents how many repetions each condition
                    end
                    resp_mat_std(k, j, i) = std(spike_rates(select));     % calculate std between trials for later DSI usage
                    %                 resp_mat_ste(k, j, i) = resp_mat_std(k, j, i)/ sqrt(length(find( (azimuth==unique_azimuth(i)) & (elevation==unique_elevation(j))&(condition_num==1) )) );
                else
                    %                resp_mat_trial{k}(t, j, i) = 0;
                    resp_mat(k, j, i) = resp_mat(k,j,1);
                    resp_mat_vector(k,j,1) =0; % for vector sum only
                    resp_mat_std(k, j, i) = 0;
                    %                 resp_mat_ste(k, j, i) = 0;
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
    %
    % calculate maximum and minimum firing rate
    max_res = max(max(max(resp_mat)));
    % max_vi_ve = max(max(max(resp_mat(2:3,:,:))));
    min_res = min(min(min(resp_mat_tran)));
    
    vector_num = length(unique_azimuth) * (length(unique_elevation)-2) + 2;
    repeat = floor( length(spike_rates) / vector_num );
    
    % % Define figure
    % xoffset=0;
    % yoffset=0;
    % figure(2);
    % set(2,'Position', [5,15 980,650], 'Name', '3D Direction Tuning');
    % orient landscape;
    % %set(0, DefaultAxesXTickMode, 'manual', DefaultAxesYTickMode, 'manual', 'DefaultAxesZTickMode', 'manual');
    % axis off;
    % % for cosine plot
    % %---------YOng's Cosine Plot------------Now disable by Katsu, 05/22/06
    % azi_cos = [1,2,3,4,5,6,7,8,9];
    % ele_sin = [-1,-0.707,0,0.707,1];
    %
    % for k=1:3 
    %     
    %     if( xoffset > 0.5)          % now temperarily 2 pictures one row and 2 one column
    %         yoffset = yoffset-0.4;
    %         xoffset = 0;
    %     end
    %     axes('position',[0.11+xoffset 0.54+yoffset 0.32 0.24]);
    % %---------Yong's Cosine Plot
    %     contourf( azi_cos, ele_sin, squeeze( resp_mat_tran(k,:,:)) );
    %     % set the same scale for visual and combined conditions but here assuming vestibular response is always smaller than that in visual and
    %     % combined conditions
    %     colorbar;
    %     % make 0 correspond to rightward and 180 correspond to leftward
    %     set(gca, 'ydir' , 'reverse');
    %     set(gca, 'xtick', [] );  
    %     set(gca, 'ytick', [] );  
    %     title( h_title{k} );
    %     
    %         % plot 1-D for mean respond as a function of elevation
    %     % notice that elevation scale is transformed by consine
    %     axes('position',[0.06+xoffset 0.54+yoffset 0.04 0.24]);
    %     for i=1:length(unique_elevation)
    %         y_elevation_mean(1,i)=mean(resp_mat_tran(k,:,i));
    %         if k ==1
    %             y_elevation_std(1,i) =std( dft_velocity([find( (elevation1==unique_elevation(i)) )]) );
    %             y_elevation_ste(1,i) =y_elevation_std(1,i)/ sqrt(length(find( (elevation1==unique_elevation(i)) )) );
    %         elseif k ==2
    %             y_elevation_std(1,i) =std( dft_acceleration([find( (elevation1==unique_elevation(i)) )]) );
    %             y_elevation_ste(1,i) =y_elevation_std(1,i)/ sqrt(length(find( (elevation1==unique_elevation(i)) )) );
    %         elseif k == 3
    %             y_elevation_std(1,i) =std( spike_rates([find( (elevation==unique_elevation(i))&(condition_num==1) )]) );
    %             y_elevation_ste(1,i) =y_elevation_std(1,i)/ sqrt(length(find( (elevation==unique_elevation(i))&(condition_num==1) )) );  
    %         end  
    %     end
    %  %---------Yong's Cosine Plot------------Now disable by Katsu, 05/22/06
    %     x_elevation=[-1,-0.707,0,0.707,1];
    %     errorbar(x_elevation,y_elevation_mean,y_elevation_ste,'ko-');%-----------Temporaly disable
    % %     errorbar(x_elevation,y_elevation_mean,y_elevation_ste,'k-');% Katsu for paper settle for m3c294
    %     xlabel('Elevation');
    %     view(90,90);
    %     set(gca, 'xtick',x_elevation);
    %     set(gca, 'ytick', []);
    %     xlim([-1, 1]);
    %     %---------Yong's Cosine Plot------------Now disable by Katsu, 05/22/06
    %     ylim([min(y_elevation_mean(1,:))-max(y_elevation_ste(1,:)), max(y_elevation_mean(1,:))+max(y_elevation_ste(1,:))]);%axis off %----------Now add axis off
    %     set(gca, 'xticklabel','-90|-45|0|45|90');
    % 
    %     
    %     % plot 1-D for mean respond as a function of azimuth
    %     axes('position',[0.11+xoffset 0.46+yoffset 0.274 0.06]);
    %     for i=1:(length(unique_azimuth) )
    %         y_azimuth_mean(1,i)=mean(resp_mat_tran(k,:,i));
    %         if k ==1
    %             y_azimuth_std(1,i) =std( dft_velocity([find( (azimuth1==unique_azimuth(i)) )]) );
    %             y_azimuth_ste(1,i) =y_azimuth_std(1,i)/ sqrt(length(find( (azimuth1==unique_azimuth(i)) )) );
    %         elseif k ==2
    %             y_azimuth_std(1,i) =std( dft_acceleration([find( (azimuth1==unique_azimuth(i)) )]) );
    %             y_azimuth_ste(1,i) =y_azimuth_std(1,i)/ sqrt(length(find( (azimuth1==unique_azimuth(i)) )) );
    %         elseif k == 3
    %             y_azimuth_std(1,i) =std( spike_rates([find( (azimuth==unique_azimuth(i))&(condition_num==1) )]) );
    %             y_azimuth_ste(1,i) =y_azimuth_std(1,i)/ sqrt(length(find( (azimuth==unique_azimuth(i))&(condition_num==1) )) );  
    %         end  
    %     end
    %     y_azimuth_mean(1,9) = mean(resp_mat_tran(k,:,1));
    %     for i=1:( length(unique_azimuth)+1 )
    %         if (i < 8)        
    %             y_azimuth_ste_tran(1,i) = y_azimuth_ste(1,8-i);
    %         elseif (i == 8)
    %             y_azimuth_ste_tran(1,i) = y_azimuth_ste(1,8);
    %         else
    %             y_azimuth_ste_tran(1,i) = y_azimuth_ste(1,7);
    %         end
    %     end
    %     x_azimuth=1:(length(unique_azimuth)+1);
    %     errorbar(x_azimuth,y_azimuth_mean,y_azimuth_ste_tran,'ko-');%----------------temporaly disable
    % %     errorbar(x_azimuth,y_azimuth_mean,y_azimuth_ste_tran,'k-');% Katsu for paper settle for m3c294
    % %     xlim( [1, length(unique_azimuth)+1] );
    %     xlim( [0.9, length(unique_azimuth)+1.1] );
    %     set(gca, 'XTickMode','manual');
    %     set(gca, 'xtick',[1,2,3,4,5,6,7,8,9]);
    %     set(gca, 'xticklabel','270|225|180|135|90|45|0|-45|-90'); %Katsu
    % %     set(gca, 'xticklabel','0|45|90|135|180|225|270|315|360'); 
    %     xlabel('Azimuth');
    %     ylim([min(y_azimuth_mean(1,:))-max(y_azimuth_ste(1,:)), max(y_azimuth_mean(1,:))+max(y_azimuth_ste(1,:))]);%axis off %----------Now add axis off
    %     xoffset=xoffset+0.48;
    % end
    
    
    %----------------------------------------------------------------------------
    % % Now show vectorsum, DSI, p and spontaneous at the top of figure
    % axes('position',[0.05,0.85, 0.9,0.1] );
    % xlim( [0,100] );
    % ylim( [0,3] );
    % h_spon = num2str(spon_resp);
    % text(0, 3, FILE);
    % %text(15,3,'Spon            Minimum        Maximum       Azi             Ele                Amp           Std             HTI             HTIerr             p');
    % text(10,3,'Protocol         p-ANOVA    p-DDI-vel    p-DDI-acc    Azi-mfr    Ele-mfr    Azi-vel    Ele-vel    Azi-acc    Ele-acc');
    % 
    % h_text= [ '     ' num2str(P_anova(1), 4) '              ' num2str(p_ddi_vel, 4) '          ' num2str(p_ddi_acc, 4) '          ' num2str(pref_az_mfr, 4) '     ' num2str(pref_el_mfr, 4) '      ' num2str(pref_az_vel, 4) '      ' num2str(pref_el_vel, 4) '      ' num2str(pref_az_acc, 4) '      ' num2str(pref_el_acc, 4)];
    % text(10,3-1,'Translation');
    % text(20,3-1, h_text );
    % axis off;
    % 
    % print -dwinc
    % close
    
    % sprint_txt = ['%s\t        %f\t'];
    % for j = 1:length(mean_rates)
    %     buff = sprintf(sprint_txt, FILE, mean_rates(j));  
    %     outfile = ['Z:\Users\Anuk\fit_analysis\data\newdata\translation_firingrates.dat'];
    %     fid = fopen(outfile, 'a');
    %     fprintf(fid, '%s', buff);
    %     fprintf(fid, '\r\n');
    %     fclose(fid);
    % end
    
    
    
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
    
    if P_anova<.05 & P_anova_vel<.05
        [r,p] = corrcoef( vel_data, mfr_data ); % vel
        vel_mfr = r(1,2);
        p_vel_mfr = p(1,2);
    end
    
    if P_anova<.05 & P_anova_acc<.05
        [r,p] = corrcoef( acc_data, mfr_data ); % vel
        acc_mfr = r(1,2);
        p_acc_mfr = p(1,2);
    end
    
    if P_anova_vel<.05 & P_anova_acc<.05
        [r,p] = corrcoef( acc_data, vel_data ); % vel
        acc_vel = r(1,2);
        p_acc_vel = p(1,2);
    end
    
    sprint_txt = ['%s\t        %f\t        %f\t       %f\t      %f\t        %f\t        %f\t'];
    buff = sprintf(sprint_txt, FILE, vel_mfr, p_vel_mfr, acc_mfr, p_acc_mfr, acc_vel, p_acc_vel);  
    outfile = ['Z:\Users\Anuk\fit_analysis\data\newdata\contour_correlation_trans_fix.dat'];
    fid = fopen(outfile, 'a');
    fprintf(fid, '%s', buff);
    fprintf(fid, '\r\n');
    fclose(fid);
    
end