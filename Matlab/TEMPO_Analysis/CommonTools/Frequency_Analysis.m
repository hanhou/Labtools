%----------------------------------------------------------------------------------------------------------------------
%Frequency_Analysis.m calculates the DFTR, phase, latency, velocity and acceleration
%components, preferred directions (for MFR, velocity and acceleration),
%DDI's, vector sum amplitudes, correlation of 3d tuning and anova for neuronal responses


% clear all;
% 
% FILE = 'm16c210r2.htb';
% 
% load(['Z:\Users\Tanya\Nodulus_frequency_analysis\all_data\' FILE '.mat']);

function Frequency_Analysis(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);

% added4RecoveringSpikeData;

disp(' '); time = clock; disp(['Start Time = ' num2str([time(4) time(5) time(6)])]);

Path_Defs;
ProtocolDefs;  %contains protocol specific keywords - 1/4/01 BJP


plot_psth = 0;
plot_spatial = 0;

dftcutoff = 2.5;    % Select the DFTR cut-off 
polswitch = 1;      % Select 1 for switching signs on the component phases
use_max_phase = 0;

% load('Z:\LabTools\Matlab\TEMPO_Analysis\ProtocolSpecific\MOOG\time_vel_acc.mat');
% load('Z:\Data\Tempo\Batch Files\Suhrud\time_vel_acc.mat');
% load('C:\MATLAB6p5\work\time_vel_acc.mat');
load('Z:\Users\Chris2\time_vel_acc.mat');

% load('Z:\LabTools\Matlab\TEMPO_Analysis\ProtocolSpecific\MOOG\bpr_lookup_table.mat');
% load('Z:\Data\Tempo\Batch Files\Suhrud\bpr_lookup_table.mat');

if Protocol == 100 || Protocol == 104
	%get the column of values for azimuth and elevation and stim_type
	temp_azimuth = data.moog_params(AZIMUTH,:,MOOG);
	temp_elevation = data.moog_params(ELEVATION,:,MOOG);
	temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG);
	temp_spike_data = data.spike_data(SpikeChan,:);
	% fixation data added by CRF, 3/2007
	temp_fix_x    =  data.moog_params(FIX_X,:,MOOG);
	temp_fix_y    =  data.moog_params(FIX_Y,:,MOOG);
	temp_fix_x(isnan(temp_fix_x)) = 0;
	temp_fix_y(isnan(temp_fix_y)) = 0;
elseif Protocol == 112
	%get the column of values for azimuth and elevation and stim_type
	temp_azimuth = data.moog_params(ROT_AZIMUTH,:,MOOG);
	temp_elevation = data.moog_params(ROT_ELEVATION,:,MOOG);
	temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG);
	temp_spike_data = data.spike_data(SpikeChan,:);
	% fixation data added by CRF, 3/2007
	temp_fix_x    =  data.moog_params(FIX_X,:,MOOG);
	temp_fix_y    =  data.moog_params(FIX_Y,:,MOOG);
	temp_fix_x(isnan(temp_fix_x)) = 0;
	temp_fix_y(isnan(temp_fix_y)) = 0;    
end
    
    
%mean firing rate of each trial depending on the start and stop offsets
temp_spike_rates = data.spike_rates(SpikeChan, :);

%--------------------------------------------------------------
% Some cells have more than the usual 26 directions, so we need to remove them
extra_directions =[];
if length(munique(temp_azimuth')) > 9  % it's 9 because the null trials (-99) are included
    unique_azimuth_withnull = [0;45;90;135;180;225;270;315;-99];%changed 9999 to 99 everywhere
    unique_elevation_withnull = [-90;-45;0;45;90;-99];
    extra_azis = ones(1,length(temp_azimuth));
    extra_eles = ones(1,length(temp_elevation));
    for n = 1:length(unique_azimuth_withnull)
        extra_azis(find(temp_azimuth == unique_azimuth_withnull(n))) = 0;
    end
    for m = 1:length(unique_elevation_withnull)
        extra_eles(find(temp_elevation == unique_elevation_withnull(m))) = 0;
    end
    extra_directions = (extra_azis | extra_eles);
    temp_elevation(extra_directions) = [];
    temp_stim_type(extra_directions) = [];
    temp_fix_x(extra_directions) = [];
    temp_fix_y(extra_directions) = [];
    temp_spike_rates(extra_directions) = [];    
    temp_azimuth(extra_directions) = [];
    
    EndTrial = EndTrial - sum(extra_directions(1:EndTrial-(BegTrial-1)));
end
%--------------------------------------------------------------

%get indices of any NULL conditions (for measuring spontaneous activity)
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

fix_x_withnull = temp_fix_x(select_trials);
fix_y_withnull = temp_fix_y(select_trials);
spike_rates_withnull = temp_spike_rates(select_trials);

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

% --------------------------------------------
% select bin size -- 50 ms (20 Hz) is default

% timebin = 50;
% time = time50; vel = vel50; acc = acc50;

timebin = 20;
time = time20; vel = vel20; acc = acc20;

% timebin = 16.6667;
% time = time16; vel = vel16; acc = acc16;

% timebin = 1;
% time = time1; vel = vel1; acc = acc1;
% --------------------------------------------

% sample frequency depends on test duration
% frequency=length(temp_spike_data)/length(select_trials); % actually its 5000 now for 5 sec 
frequency = 5000;

% length of x-axis
x_length = round(frequency/timebin);
% x-axis for plot PSTH
x_time=1:x_length;
x_time_bincenter = timebin/2000;            % using bin-centers now: start at time = 0, plus a half-bin...
for nn = 1 : round(5/(timebin/1000)) - 1    % then fill it out with a full 2-seconds (i.e., 40 bins with timebin = 50 ms)
    x_time_bincenter(end+1) = x_time_bincenter(end) + (timebin/1000);
end

% remove 'extra direction' trials from spikedata
discard_trials = find(extra_directions);
for i = 1 : length(discard_trials)
    temp_spike_data( 1, ((discard_trials(i)-1)*frequency +1) : discard_trials(i)*frequency ) = 99;
end
temp_spike_data = temp_spike_data( temp_spike_data~=99 );

% remove trials outside Begtrial~Endtrial
discard_trials = [find(trials <BegTrial | trials >EndTrial)];
for i = 1 : length(discard_trials)
    temp_spike_data( 1, ((discard_trials(i)-1)*frequency +1) : discard_trials(i)*frequency ) = 99;
end
spike_data_withnull = temp_spike_data( temp_spike_data~=99 );

% remove null trials and trials outside Begtrial~Endtrial
discard_trials = [find(null_trials==1 | trials <BegTrial | trials >EndTrial)];
for i = 1 : length(discard_trials)
    temp_spike_data( 1, ((discard_trials(i)-1)*frequency +1) :  discard_trials(i)*frequency ) = 99;
end
spike_data = temp_spike_data( temp_spike_data~=99 );

% count spikes from raster data (spike_data)
max_count = 1;
time_step = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
current_stim = 2;  % vestibular=1, visual=2
actual_stim = current_stim;

% For visual-only blocks, the index 'k' (below) will only reach 1, and no data
% will be created at index k = 2.  Therefore must change current_stim to 1. -CRF
if unique_condition_num == 2 & current_stim == 2
    current_stim = 1;
    actual_stim = 2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k=1:length(unique_condition_num)  % it's the stim_type
    for j=1:length(unique_elevation)
        for i=1:length(unique_azimuth)         % SELECT now restricted to only 0 deg fixation trials -- CRF 3/2007
            select = logical(azimuth==unique_azimuth(i) & elevation==unique_elevation(j) & condition_num==unique_condition_num(k) & fix_x==0 & fix_y==0 );            
            if sum(select) > 0
                resp{k}(j,i) = mean(spike_rates(select));
                act_found = find( select==1 );
                % count spikes per timebin on every same condition trials
                clear temp_count dummy_count;
                for repeat=1:length(act_found) 
                    for n=1:(x_length)
                        temp_count(repeat,n)=sum( spike_data( 1, frequency*(act_found(repeat)-1) + time_step : round(frequency*(act_found(repeat)-1)+n*timebin) ) );
                        dummy_count{repeat}(n) = sum(spike_data(1,(frequency*(act_found(repeat)-1)+time_step):round(frequency*(act_found(repeat)-1)+n*timebin)));
                        if timebin == 16.6667
                            if floor((n+1)/3) == (n+1)/3
                                time_step = floor(time_step + timebin);
						    else                               % kluge for 16.6667 ms timebin (just rounding was missing bins in spike_data)
                                time_step = round(time_step + timebin);
							end
                        else
                            time_step = time_step + timebin;   % this may need to be modified for timebins other than 50
                        end
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

% Do the same, but for null trials (for new significance test on DFTR below) -- CRF 10/2007
select = null_trials(select_trials) & temp_fix_x(select_trials)==0 & temp_fix_y(select_trials)==0;
if sum(select) > 0
    resp_null = mean(spike_rates_withnull(select));
    act_found = find( select==1 );
    % count spikes per timebin on every same condition trials
    clear temp_count dummy_count;
    for repeat=1:length(act_found) 
        for n=1:(x_length)
            temp_count(repeat,n)=sum( spike_data_withnull( 1, frequency*(act_found(repeat)-1) + time_step : round(frequency*(act_found(repeat)-1)+n*timebin) ) );
            dummy_count{repeat}(n) = sum(spike_data_withnull(1,(frequency*(act_found(repeat)-1)+time_step):round(frequency*(act_found(repeat)-1)+n*timebin)));
            if timebin == 16.6667
                if floor((n+1)/3) == (n+1)/3
                    time_step = floor(time_step + timebin);
			    else                               % kluge for 16.6667 ms timebin (just rounding was missing bins in spike_data)
                    time_step = round(time_step + timebin);
				end
            else
                time_step = time_step + timebin;   % this may need to be modified for timebins other than 50
            end
        end
        time_step=1;
    end
    count_y_trial_null = temp_count;
    dim=size(temp_count);
    if dim(1) > 1;
        count_y_null = mean(temp_count);
    else
        count_y_null = temp_count;     % for only one repetition cases
    end
else
    resp_null=0; 
    count_y_null=0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% order in which the directions are plotted
plot_col = [1 1 1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 6 6 6 7 7 7 8 8 8];
plot_row = [5 4 3 2 1 4 3 2 4 3 2 4 3 2 4 3 2 4 3 2 4 3 2 4 3 2];

% find when in the 5-sec trial period that the spike data cuts off
sum_dat = zeros(1,x_length);
for j = 1:26
    sum_dat = sum_dat + count_y{plot_col(j), plot_row(j), current_stim};
end
% and accordingly set length of data to read from PSTH
for ii = 1:x_length
    if sum(sum_dat(ii:end)) == 0
        span = ii;
      break
  else
       span = length(count_y{plot_col(j), plot_row(j), current_stim});% this line was added at 06/16/10 by Chris reccomendation
    end
end %before two ends were after breat
        
if span < .80*x_length
    span = round(.68*x_length); % usually this is how much data there is
else % this is one of katsu's cells
    span = round(.80*x_length);
end

% ----------------------------------------------------------------------------------------------
% Define some time points for segmenting the response into the segment of interest -CRF 6/2007
% ----------------------------------------------------------------------------------------------

% The stimulus velocity peaks at time t = 2092 ms, however the 'onset' is actually t = 996 ms, according to the event code* (4).
% Nevertheless, for simplicity we will assume a 2-second stimulus window beginning at t = 1092 ms and ending at t = 3092.
% *(Oddly, the stim offset code (5) occurs at t = 3006, so there's an extra 10 ms in there somehow.)

% Take baseline data as the 200 ms before stimulus onset (t = 996 ms)
tempbins = find(x_time_bincenter > .796);
baseline_begin = tempbins(1);
tempbins = find(x_time_bincenter < .996);
baseline_end = tempbins(end);

% Data for hilbert transform (temporal envelope) can be the full 2-second window
tempbins = find(x_time_bincenter > 1.092);
hilbert_begin = tempbins(1);
tempbins = find(x_time_bincenter < 3.092);
hilbert_end = tempbins(end);

% 
% 
% % TEMP: expanded window to address 'edge effect' concern
% hilbert_begin = hilbert_begin - (span-hilbert_end);
% hilbert_end = span;
% 
% 

% Data from 1400 to 2700 ms into stim period, from which the peak is taken
% (but only as a failsafe for when center of mass doesn't fall in this range, which is relatively rare)
tempbins = find(x_time_bincenter > 1.492);
peak_begin = tempbins(1);
tempbins = find(x_time_bincenter < 2.792);
peak_end = tempbins(end);

% Finally, find the bin closest to the stimulus peak velocity (2.092 s)
peak_stim = find(abs(x_time_bincenter - 2.092) == min(abs(x_time_bincenter - 2.092)));
peak_stim = peak_stim(1);

for j = 1:27
    emptyerror1(j) = 0;
    emptyerror2(j) = 0;
end


% -----------------------------------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------------------------------
% perform frequency analysis---------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------------------------------
for j=1:27
    if j == 27
        gdat = count_y_null;
    else
        gdat = count_y{plot_col(j), plot_row(j), current_stim};
    end
   
    %     NEW: baseline diff post-pre (see single trials loop below) to plot spatial tuning of position signal
    pre_baseline_mean_avg(j) = sum(gdat(baseline_begin:baseline_end)) / (timebin/1000);
    post_baseline_mean_avg(j) = sum(gdat(hilbert_end+5:hilbert_end+14)) / (timebin/1000);  % hard coded for now
    baseline_diff_avg(j) = post_baseline_mean_avg(j) - pre_baseline_mean_avg(j);
        
    % %  %%NEW STUFF APR 15th BEGIN  
    datatogetmean = gdat;
    meanfiring = mean(datatogetmean(30:60));
    
    sizeofbins = 10;
    slider = 1;
    startat = hilbert_begin;
    nn = 1;
    endat = startat +sizeofbins-1;
    while endat < hilbert_end
        if sum(datatogetmean(startat : endat) )~=0 
            meandata(nn) = mean(datatogetmean(startat : endat) );
        else 
            meandata(nn) = 0;
        end
        startat = startat + slider;
        endat = startat + sizeofbins -1;
        nn = nn + 1;
    end
    
    maxmean = meandata(find(meandata == max(meandata)));
    bpr(j) = meanfiring/maxmean(1);
    
% % %     figure(j+1); bar(x_time, gdat);
% % %     hold on; plot(x_time(slider*find(meandata == max(meandata)))+hilbert_begin, maxmean, 'r*');
% % end
% %     
% % %%NEW STUFF APR 15 END
% %     
    total_data{j} = gdat(1:span);
    %---------------------------------------
    % find hilbert transform to get peak.
    shift_win(j) = 0;
    flagged_dir(j) = 0;
    % running mean
    binsize = 3; % take mean of 3 neighbour bins
    slidesize = 1; % with a slide of 1 bin
    nn = 1;
    start_bin = 1;
    end_bin = start_bin + binsize -1;
    while end_bin < x_length % the length of gdat
        if sum( gdat(start_bin : end_bin) )~=0 % avoid divided by 0 later
            runmean(nn) = mean( gdat(start_bin : end_bin) );
        else 
            runmean(nn) = 0;
        end
        start_bin = start_bin + slidesize;
        end_bin = start_bin + binsize -1;
        nn = nn + 1;
    end
    hildata = zeros(1,x_length);
    hildata(2:length(runmean)+1) = runmean;
    dummy_dat{j} = hildata(hilbert_begin : hilbert_end) - mean(hildata(baseline_begin : baseline_end));
    x_time_bincenter_hil = x_time_bincenter(hilbert_begin : hilbert_end);
    hil{j} = hilbert(dummy_dat{j});
    abb = abs(hil{j});
%     HT = abb - mean(abb(round(.01*x_length) : round(.03*x_length))); % remove the mean before finding the center of mass or else it will be biased (???)
    HT = abb - min(abb);
    
    % Find four largest bins as candidates for the peak of the envelope.
    % The bin with the largest mean over 3 bins will be the peak.
    % This is only used when the center of mass falls outside the range 1.5s to 2.7s (not often)
    dum_HT = HT(peak_begin-hilbert_begin+1 : peak_end-hilbert_begin+1);
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
    dum_HT = HT(peak_begin-hilbert_begin+1 : peak_end-hilbert_begin+1);
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
    peak_temp = m(maxbin) + peak_begin - 1;

    H_time = x_time_bincenter(hilbert_begin : hilbert_end);
    CMsum = HT .* H_time;
    if sum(HT) > 0
        center_mass = sum(CMsum)/sum(HT);
        peak_t = find( abs(x_time_bincenter - center_mass) == min(abs(x_time_bincenter - center_mass)) );% point which is closest to the center of mass is used.
        if x_time_bincenter(peak_t) >= 1.4 & x_time_bincenter(peak_t) <= 2.7
            peak_time(j) = peak_t(1);
        else 
            peak_time(j) = peak_temp(1);
            flagged_dir(j) = 1;
        end
    else
        peak_time(j) = peak_stim;
    end

    shift_win(j) = peak_time(j) - peak_stim; % 'delay' (shift of window for fourier transform to get rid of phase shift caused by delayed response)
%     shift_win(j) = 10;  % temporarily fix delay
    
    %if the window is being shifted to extend into the region where the
    %recording has stopped, then we will randomly sample from the front of
    %the response and fill out the window.
    signal = total_data{j}(baseline_begin : baseline_end);
    span_new(j) = span;
    
    if hilbert_end + shift_win(j) > span
        for i = 1:(hilbert_end + shift_win(j) - span)
            rand_sig = randperm(length(signal));
            total_data{j}(span+i) = signal(rand_sig(1));
            span_new(j) = span+i;
        end
    end
    
    % window is shifted depending on the delay of the peak of the envelope (from the hilbert trans), relative to the peak of the stimulus
    gauss_dat{j} = total_data{j}(hilbert_begin+shift_win(j) : hilbert_end+shift_win(j));
%     gauss_dat{j} = total_data{j}(hilbert_begin+shift_win(j)-1 : hilbert_end+shift_win(j)-1);      %******************************************
    gauss_time{j} = x_time_bincenter(hilbert_begin+shift_win(j) : hilbert_end+shift_win(j));        % NOT SURE HERE, BECAUSE SEGMENT HAS 40 BINS, AND
%     gauss_time{j} = x_time_bincenter(hilbert_begin+shift_win(j)-1 : hilbert_end+shift_win(j)-1);  % THUS CANNOT BE CENTERED ON ONE PARTICULAR BIN.
                                                                                                    % FOR 50-MS TIMEBIN, THIS CAUSES A PHASE DIFF OF 9 DEG.
                                                                                                    % ******************************************
    % Also recompute MFR with the delay-shifted window
    if j == 27
        resp2_null = sum(total_data{j}(hilbert_begin+shift_win(j) : hilbert_end+shift_win(j))) / 2;
    else
        resp2{current_stim}(plot_row(j),plot_col(j)) = sum(total_data{j}(hilbert_begin+shift_win(j) : hilbert_end+shift_win(j))) / 2;
    end
                                                                                                    
%     time_short = time(hilbert_begin:hilbert_end);
%     vel_short = vel(hilbert_begin:hilbert_end);
%     acc_short = acc(hilbert_begin:hilbert_end);

    %calculate DFT ratio
    %40 pt DFT with mean removed
    [f, amp, resp_phase] = FT(gauss_time{j}, gauss_dat{j}, length(gauss_time{j}), 1, 0);
    f = round(f*100)/100; % get rid of some floating point issues
    
    % NOTE DFT cut-off added. CRF and SR, 3-18-2008.
    
    f1 = mean(amp(find(f > 0 & f <= dftcutoff)));
    f2 = mean(amp(find(f > dftcutoff)));
    if f2 == 0
        fourier_ratio(j) = 0;
    else        %DFTR
        fourier_ratio(j) = f1/f2;
    end

    %phase of response
    max_c = find(amp == max(amp(find(f > 0 & f <= dftcutoff))));
    max_comp(j) = max_c(1);
    max_p = mean(resp_phase(max_c));
    max_phase(j) = max_p(1);

    %  New approach to vel and accel DFTR: multiply the amplitude of each frequency component by the 
	%  cosine (vel) or sine (acc) of its phase angle (essentially scaling each by a measure of how
	%  close its phase is to the phase of the vel and accel stimulus profiles), then recompute a separate 
	%  vel-DFTR and acc-DFTR from these scaled amplitudes.     GCD and CRF, 7-20-2007

    if (use_max_phase == 0)
        pol = -1;
        for s = find(f > 0 & f <= dftcutoff)
            amp_vel(s) = pol*amp(s)*cos(resp_phase(s));
            amp_acc(s) = pol*amp(s)*sin(resp_phase(s));
            if polswitch
                pol = -pol;
            end
        end
        v1 = mean(amp_vel(find(f > 0 & f <= dftcutoff)));
        if f2 == 0
            fourier_ratio_vel(j) = 0;
        else
            fourier_ratio_vel(j) = v1/f2;
        end
        a1 = mean(amp_acc(find(f > 0 & f <= dftcutoff)));
        if f2 == 0
            fourier_ratio_acc(j) = 0;
        else
            fourier_ratio_acc(j) = a1/f2;
        end
    elseif (use_max_phase == 1)
        if polswitch
	        if max_comp(j) == 3
	            pol = 1;
            else
	            pol = -1;
	        end
        else
            pol = -1;
        end
        fourier_ratio_vel(j) = pol*cos(max_phase(j))*f1 / f2;
        fourier_ratio_acc(j) = pol*sin(max_phase(j))*f1 / f2;
    end
    
    % Define a vel and acc 'percentage'
    dft_sum(j) = abs(fourier_ratio_vel(j)) + abs(fourier_ratio_acc(j));
    vel_pct(j) = fourier_ratio_vel(j) / dft_sum(j);
    vel_pct(j) = round(vel_pct(j)*100) / 100;
    acc_pct(j) = fourier_ratio_acc(j) / dft_sum(j);
    acc_pct(j) = round(acc_pct(j)*100) / 100;

%     % Adjust vel/acc percentages using a lookup-table approach, based on fake data
%     bpr(j) = round(bpr(j)*100) / 100;
%     vel_pct_adj_set = [];
%     acc_pct_adj_set = [];
%     emptyerror1(j) = 0;
%     emptyerror2(j) = 0;
%     
%     if isempty(find(bpr_lookup == bpr(j) & vel_pct_output == vel_pct(j)))
%         emptyerror1(j) = 1;
%         vel_pct_adj_set = vel_pct(j);
%     else
%         vel_pct_adj_set = vel_pct_input( find(bpr_lookup == bpr(j) & vel_pct_output == vel_pct(j)) );
% 	end
% 	
%     if isempty(find(bpr_lookup == bpr(j) & acc_pct_output == acc_pct(j)))
%         emptyerror2(j) = 1;
%         acc_pct_adj_set = acc_pct(j);
%     else
%         acc_pct_adj_set = acc_pct_input( find(bpr_lookup == bpr(j) & acc_pct_output == acc_pct(j)) );
% 	end
%     
%     index = 1;
%     for a = 1:length(vel_pct_adj_set)
%         for b = 1:length(acc_pct_adj_set)
%             if abs(vel_pct_adj_set(a)) + abs(acc_pct_adj_set(b)) == 1
%                 vel_pct_good(index) = vel_pct_adj_set(a);
%                 acc_pct_good(index) = acc_pct_adj_set(b);
%                 index = index + 1;
%             end
%         end
%     end
% 
%     if emptyerror1(j)
%         vel_pct_adj(j) = vel_pct(j);
%     else
%         vel_pct_adj(j) = mean(unique(vel_pct_good));
%     end
%     if emptyerror2(j)
%         acc_pct_adj(j) = acc_pct(j);
%     else
%         acc_pct_adj(j) = mean(unique(acc_pct_good));
%     end
%     
%     % Lastly, recompute vel and acc DFTRs via the new _pcts
%     vel_dftr_adj(j) = vel_pct_adj(j) * dft_sum(j);
%     acc_dftr_adj(j) = acc_pct_adj(j) * dft_sum(j);


% % % don't adjust dftr's for now
% %     
% %     fourier_ratio_vel(j) = vel_dftr_adj(j);
% %     fourier_ratio_acc(j) = acc_dftr_adj(j);
% % 
    % instead just assign them as non-adjusted
    vel_dftr_adj(j) = fourier_ratio_vel(j);
    acc_dftr_adj(j) = fourier_ratio_acc(j);
    vel_pct_adj(j) = vel_pct(j);
    acc_pct_adj(j) = acc_pct(j);
    
    %calculate p-values
%     [dft_p(j) dft_p_vel(j) dft_p_acc(j)] = DFT_perm(gauss_time{j}, gauss_dat{j}, fourier_ratio(j), fourier_ratio_vel(j), fourier_ratio_acc(j), 1000, dftcutoff);
      dft_p(j) = 999; dft_p_vel(j) = 999; dft_p_acc(j) = 999;

end
% end frequency analysis-------------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------------------------------

% new test: bootstrap single trials to test DFTR of each dir vs. null-trial DFTR -- CRF 3-26-10
num_boots = 1000
for j = 1:26
    j
    fourier_ratio_boot = zeros(1,num_boots);
    if fourier_ratio(j) > fourier_ratio(27) % only need to do it if the raw value is > Null
        for b = 1:num_boots
            for i = 1:length(count_y_trial{plot_col(j), plot_row(j), current_stim}(:,1)) 
                randboot = randperm(length(count_y_trial{plot_col(j), plot_row(j), current_stim}(:,1)));
                count_y_boot(i,:) = count_y_trial{plot_col(j), plot_row(j), current_stim}(randboot(1),:);
            end
            
            gdat = mean(count_y_boot);  % ************** NOW RUN FA: **************
            
            %     NEW: baseline diff post-pre (see single trials loop below) to plot spatial tuning of position signal
            pre_baseline_mean_avg_boot(j) = sum(gdat(baseline_begin:baseline_end)) / (timebin/1000);
            post_baseline_mean_avg_boot(j) = sum(gdat(hilbert_end+5:hilbert_end+14)) / (timebin/1000);  % hard coded for now
            baseline_diff_avg_boot(j) = post_baseline_mean_avg_boot(j) - pre_baseline_mean_avg_boot(j);
                
            % %  %%NEW STUFF APR 15th BEGIN  
            datatogetmean = gdat;
            meanfiring = mean(datatogetmean(30:60));
            
            sizeofbins = 10;
            slider = 1;
            startat = hilbert_begin;
            nn = 1;
            endat = startat +sizeofbins-1;
            while endat < hilbert_end
                if sum(datatogetmean(startat : endat) )~=0 
                    meandata(nn) = mean(datatogetmean(startat : endat) );
                else 
                    meandata(nn) = 0;
                end
                startat = startat + slider;
                endat = startat + sizeofbins -1;
                nn = nn + 1;
            end
            
            maxmean = meandata(find(meandata == max(meandata)));
            bpr_boot(j) = meanfiring/maxmean(1);
            
		% % %     figure(j+1); bar(x_time, gdat);
		% % %     hold on; plot(x_time(slider*find(meandata == max(meandata)))+hilbert_begin, maxmean, 'r*');
		% % end
		% %     
		% % %%NEW STUFF APR 15 END
		% %     
            total_data_boot{j} = gdat(1:span);
            %---------------------------------------
            % find hilbert transform to get peak.
            shift_win_boot(j) = 0;
            flagged_dir_boot(j) = 0;
            % running mean
            binsize = 3; % take mean of 3 neighbour bins
            slidesize = 1; % with a slide of 1 bin
            nn = 1;
            start_bin = 1;
            end_bin = start_bin + binsize -1;
            while end_bin < x_length % the length of gdat
                if sum( gdat(start_bin : end_bin) )~=0 % avoid divided by 0 later
                    runmean(nn) = mean( gdat(start_bin : end_bin) );
                else 
                    runmean(nn) = 0;
                end
                start_bin = start_bin + slidesize;
                end_bin = start_bin + binsize -1;
                nn = nn + 1;
            end
            hildata = zeros(1,x_length);
            hildata(2:length(runmean)+1) = runmean;
            dummy_dat_boot{j} = hildata(hilbert_begin : hilbert_end) - mean(hildata(baseline_begin : baseline_end));
            x_time_bincenter_hil = x_time_bincenter(hilbert_begin : hilbert_end);
            hil_boot{j} = hilbert(dummy_dat_boot{j});
            abb = abs(hil_boot{j});
		%     HT = abb - mean(abb(round(.01*x_length) : round(.03*x_length))); % remove the mean before finding the center of mass or else it will be biased (???)
            HT = abb - min(abb);
            
            % Find four largest bins as candidates for the peak of the envelope.
            % The bin with the largest mean over 3 bins will be the peak.
            % This is only used when the center of mass falls outside the range 1.5s to 2.7s (not often)
            dum_HT = HT(peak_begin-hilbert_begin+1 : peak_end-hilbert_begin+1);
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
            dum_HT = HT(peak_begin-hilbert_begin+1 : peak_end-hilbert_begin+1);
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
            peak_temp = m(maxbin) + peak_begin - 1;
		
            H_time = x_time_bincenter(hilbert_begin : hilbert_end);
            CMsum = HT .* H_time;
            if sum(HT) > 0
                center_mass = sum(CMsum)/sum(HT);
                peak_t = find( abs(x_time_bincenter - center_mass) == min(abs(x_time_bincenter - center_mass)) );% point which is closest to the center of mass is used.
                if x_time_bincenter(peak_t) >= 1.4 & x_time_bincenter(peak_t) <= 2.7
                    peak_time_boot(j) = peak_t(1);
                else 
                    peak_time_boot(j) = peak_temp(1);
                    flagged_dir_boot(j) = 1;
                end
            else
                peak_time_boot(j) = peak_stim;
            end

            shift_win_boot(j) = peak_time_boot(j) - peak_stim; % 'delay' (shift of window for fourier transform to get rid of phase shift caused by delayed response)
		%     shift_win_boot(j) = 10;  % temporarily fix delay
            
            %if the window is being shifted to extend into the region where the
            %recording has stopped, then we will randomly sample from the front of
            %the response and fill out the window.
            signal = total_data_boot{j}(baseline_begin : baseline_end);
            span_new_boot(j) = span;
            
            if hilbert_end + shift_win_boot(j) > span
                for i = 1:(hilbert_end + shift_win_boot(j) - span)
                    rand_sig = randperm(length(signal));
                    total_data_boot{j}(span+i) = signal(rand_sig(1));
                    span_new_boot(j) = span+i;
                end
            end
            
            % window is shifted depending on the delay of the peak of the envelope (from the hilbert trans), relative to the peak of the stimulus
            gauss_dat_boot{j} = total_data_boot{j}(hilbert_begin+shift_win_boot(j) : hilbert_end+shift_win_boot(j));
		%     gauss_dat_boot{j} = total_data_boot{j}(hilbert_begin+shift_win_boot(j)-1 : hilbert_end+shift_win_boot(j)-1);      %******************************************
            gauss_time_boot{j} = x_time_bincenter(hilbert_begin+shift_win_boot(j) : hilbert_end+shift_win_boot(j));        % NOT SURE HERE, BECAUSE SEGMENT HAS 40 BINS, AND
		%     gauss_time_boot{j} = x_time_bincenter(hilbert_begin+shift_win_boot(j)-1 : hilbert_end+shift_win_boot(j)-1);  % THUS CANNOT BE CENTERED ON ONE PARTICULAR BIN.
                                                                                                            % FOR 50-MS TIMEBIN, THIS CAUSES A PHASE DIFF OF 9 DEG.
            %calculate DFT ratio
            %40 pt DFT with mean removed
            [f, amp, resp_phase] = FT(gauss_time_boot{j}, gauss_dat_boot{j}, length(gauss_time_boot{j}), 1, 0);
            f = round(f*100)/100; % get rid of some floating point issues
            
            f1 = mean(amp(find(f > 0 & f <= dftcutoff)));
            f2 = mean(amp(find(f > dftcutoff)));
            if f2 == 0
                fourier_ratio_boot(b) = 0;
            else        %DFTR
                fourier_ratio_boot(b) = f1/f2;
            end
            
            % ********** END FA FOR BOOTSTRAP ***********
            
        end
        dft_p_vs_Null(j) = 1 - sum(fourier_ratio_boot > fourier_ratio(27))/num_boots;
    else
        dft_p_vs_Null(j) = 99;
    end
end
dft_p_vs_Null(27) = NaN;


% % ----------------------------------------------------------------------------
% % Calculate the preferred direction using vector sum
% % ----------------------------------------------------------------------------
DFT_ratio = zeros(5,8); vel_sum = zeros(5,8); acc_sum = zeros(5,8);
for i = 1:26
    DFT_ratio(plot_row(i), plot_col(i)) = fourier_ratio(i);
    vel_sum(plot_row(i), plot_col(i)) = fourier_ratio_vel(i);
    acc_sum(plot_row(i), plot_col(i)) = fourier_ratio_acc(i);
end

[mfr_pref_az mfr_pref_el mfr_pref_amp] = vectorsum(resp2{current_stim});
[pref_az_full pref_el_full pref_amp_full] = vectorsum(DFT_ratio);  % added pref dir based on full (raw) DFTR  -CRF
[pref_az_vel pref_el_vel pref_amp_vel] = vectorsum(vel_sum);
[pref_az_acc pref_el_acc pref_amp_acc] = vectorsum(acc_sum);

%----------------------------------------------------------------------------
% make contour maps to compare with averaged single-trial versions below ('Spatial Correlation') -CRF
DFT_ratio_temp = DFT_ratio;
DFT_ratio_temp(1,:) = DFT_ratio_temp(1,1);
DFT_ratio_temp(5,:) = DFT_ratio_temp(5,1);
DFT_ratio_tran(:,1) = DFT_ratio_temp(:,7);
DFT_ratio_tran(:,2) = DFT_ratio_temp(:,6);
DFT_ratio_tran(:,3) = DFT_ratio_temp(:,5);
DFT_ratio_tran(:,4) = DFT_ratio_temp(:,4);
DFT_ratio_tran(:,5) = DFT_ratio_temp(:,3);
DFT_ratio_tran(:,6) = DFT_ratio_temp(:,2);
DFT_ratio_tran(:,7) = DFT_ratio_temp(:,1);
DFT_ratio_tran(:,8) = DFT_ratio_temp(:,8);
DFT_ratio_tran(:,9) = DFT_ratio_tran(:,1);
DFT_ratio_tran = flipud(DFT_ratio_tran);

vel_sum_temp = vel_sum;
vel_sum_temp(1,:) = vel_sum_temp(1,1);
vel_sum_temp(5,:) = vel_sum_temp(5,1);
vel_sum_tran(:,1) = vel_sum_temp(:,7);
vel_sum_tran(:,2) = vel_sum_temp(:,6);
vel_sum_tran(:,3) = vel_sum_temp(:,5);
vel_sum_tran(:,4) = vel_sum_temp(:,4);
vel_sum_tran(:,5) = vel_sum_temp(:,3);
vel_sum_tran(:,6) = vel_sum_temp(:,2);
vel_sum_tran(:,7) = vel_sum_temp(:,1);
vel_sum_tran(:,8) = vel_sum_temp(:,8);
vel_sum_tran(:,9) = vel_sum_tran(:,1);
vel_sum_tran = flipud(vel_sum_tran);

acc_sum_temp = acc_sum;
acc_sum_temp(1,:) = acc_sum_temp(1,1);
acc_sum_temp(5,:) = acc_sum_temp(5,1);
acc_sum_tran(:,1) = acc_sum_temp(:,7);
acc_sum_tran(:,2) = acc_sum_temp(:,6);
acc_sum_tran(:,3) = acc_sum_temp(:,5);
acc_sum_tran(:,4) = acc_sum_temp(:,4);
acc_sum_tran(:,5) = acc_sum_temp(:,3);
acc_sum_tran(:,6) = acc_sum_temp(:,2);
acc_sum_tran(:,7) = acc_sum_temp(:,1);
acc_sum_tran(:,8) = acc_sum_temp(:,8);
acc_sum_tran(:,9) = acc_sum_tran(:,1);
acc_sum_tran = flipud(acc_sum_tran);

mfr_sum_temp = resp2{current_stim};
mfr_sum_temp(1,:) = mfr_sum_temp(1,1);
mfr_sum_temp(5,:) = mfr_sum_temp(5,1);
mfr_sum_tran(:,1) = mfr_sum_temp(:,7);
mfr_sum_tran(:,2) = mfr_sum_temp(:,6);
mfr_sum_tran(:,3) = mfr_sum_temp(:,5);
mfr_sum_tran(:,4) = mfr_sum_temp(:,4);
mfr_sum_tran(:,5) = mfr_sum_temp(:,3);
mfr_sum_tran(:,6) = mfr_sum_temp(:,2);
mfr_sum_tran(:,7) = mfr_sum_temp(:,1);
mfr_sum_tran(:,8) = mfr_sum_temp(:,8);
mfr_sum_tran(:,9) = mfr_sum_tran(:,1);
mfr_sum_tran = flipud(mfr_sum_tran);

if plot_spatial
	figure; set(gcf, 'position', [50, 200, 1200, 350]);
    subplot(1,4,1); contourf(DFT_ratio_tran); set(gca, 'ytick', [1 2 3 4 5] ); set(gca, 'xtick', [1 2 3 4 5 6 7 8 9] ); set(gca, 'xticklabel','270|225|180|135|90|45|0|-45|-90'); set(gca, 'yticklabel','90|45|0|-45|-90'); title([FILE '  DFTR']); colorbar;
    subplot(1,4,2); contourf(vel_sum_tran); set(gca, 'ytick', [1 2 3 4 5] ); set(gca, 'xtick', [1 2 3 4 5 6 7 8 9] ); set(gca, 'xticklabel','270|225|180|135|90|45|0|-45|-90'); set(gca, 'yticklabel','90|45|0|-45|-90'); title([FILE '  DFTR Vel']); colorbar;
	subplot(1,4,3); contourf(acc_sum_tran); set(gca, 'ytick', [1 2 3 4 5] ); set(gca, 'xtick', [1 2 3 4 5 6 7 8 9] ); set(gca, 'xticklabel','270|225|180|135|90|45|0|-45|-90'); set(gca, 'yticklabel','90|45|0|-45|-90'); title([FILE '  DFTR Acc']); colorbar;
	subplot(1,4,4); contourf(mfr_sum_tran); set(gca, 'ytick', [1 2 3 4 5] ); set(gca, 'xtick', [1 2 3 4 5 6 7 8 9] ); set(gca, 'xticklabel','270|225|180|135|90|45|0|-45|-90'); set(gca, 'yticklabel','90|45|0|-45|-90'); title([FILE '  MFR']); colorbar;
	orient landscape; set(gcf,'PaperPositionMode', 'auto');
%	 print;
%	 saveas(gcf,['C:\MATLAB6p5\work\Freq_analysis_figures\' FILE '.fig']);
%	 close;
end

if plot_psth
%--------------------------------------------------------------------------
% Manually plot circular array of PSTH's (at each elev) -- CRF 6/19/07
for q = -90:45:90
	j = find(unique_elevation==q);
    k = current_stim;
	x_start = [StartEventBin(1,1)/timebin, StartEventBin(1,1)/timebin];
	x_stop =  [StopEventBin(1,1)/timebin,  StopEventBin(1,1)/timebin];
	y_marker=[0,max_count];
	figure; set(gcf,'Position', [30+(q+45)*6,50+(q+45) 700,700]);
	title([FILE ' - Elev = ' num2str(q)]); axis off;
	xscale = .28;
	yscale = .32;
	xoffset = .42;
	yoffset = .42;
	for i=1:length(unique_azimuth)
        axes('position',[xscale*cos(unique_azimuth(i)*pi/180)+xoffset yscale*sin(unique_azimuth(i)*pi/180)+yoffset 0.18 0.12]);
        bar(x_time, count_y{i,j,k});
        text(0.15, 0.9*max_count, num2str(find(plot_col==i & plot_row==j)));
        hold on;
        plot( x_start, y_marker, 'k-');
        plot( x_stop,  y_marker, 'k-');
        set( gca, 'xticklabel', '0|1|2|3' );
        set( gca, 'yticklabel', ' ' );
        % set the same scale for all plot
        xlim([0,span]);
        ylim([0,max_count]);
        % add vel and acc profiles, shifted by delay
        delay = shift_win(find(plot_col==i & plot_row==find(unique_elevation==q)));
        Vel = vel; Acc = acc;
        if delay > 0
            Vel(1+delay:end+delay) = vel;
            Vel(end+1-delay:end) = [];
            Acc(1+delay:end+delay) = acc;
            Acc(end+1-delay:end) = [];
        elseif delay < 0
            Vel(end+1:end+abs(delay)) = 0;
            Vel = Vel(1+abs(delay):end);
            Acc(end+1:end+abs(delay)) = 0;
            Acc = Acc(1+abs(delay):end);
        end
        plot(x_time,(Vel/max(Vel))*0.9*max_count,'g');
        plot(x_time,(Acc/max(Acc))*0.4*max_count+max_count/2,'r')
        vel_p = dft_p_vel(find(plot_col==i & plot_row==find(unique_elevation==q)));
        acc_p = dft_p_acc(find(plot_col==i & plot_row==find(unique_elevation==q)));
% % %             hilb1 = zeros(1,22);
% % %             hilb3 = zeros(1,38);
        hilb1 = zeros(1, hilbert_begin-1);
        hilb3=zeros(1,(x_time(end)-hilbert_end));
%         if q==-90
%             hilbert_to_plot = abs(hil{find(plot_col==1 & plot_row==find(unique_elevation==q))});  %added plot of the hilbert tranform -- SR 02/26/08
%             hilb2 = hilbert_to_plot - min(hilbert_to_plot);
%             hilbert_plot = cat(2, hilb1, hilb2, hilb3);
%             plot(x_time, hilbert_plot, 'y');
%         elseif q==90
%             hilbert_to_plot = abs(hil{find(plot_col==1 & plot_row==find(unique_elevation==q))});  %added plot of the hilbert tranform -- SR 02/26/08
%             hilb2 = hilbert_to_plot - min(hilbert_to_plot);
%             hilbert_plot = cat(2, hilb1, hilb2, hilb3);
%             plot(x_time, hilbert_plot, 'y');
%         else
%             hilbert_to_plot = abs(hil{find(plot_col==i & plot_row==find(unique_elevation==q))});  %added plot of the hilbert tranform -- SR 02/26/08
%             hilb2 = hilbert_to_plot - min(hilbert_to_plot);
%             hilbert_plot = cat(2, hilb1, hilb2, hilb3);
%             plot(x_time, hilbert_plot, 'y');
%         end
        title(['azi = ' num2str(unique_azimuth(i)) ' delay = ' num2str(delay*timebin/1000)]);
        xlabel(['v=' num2str(vel_sum(find(unique_elevation==q),i),2) ',' num2str(vel_p,2) ' a=' num2str(acc_sum(find(unique_elevation==q), i),2) ',' num2str(acc_p,2) 'df=' num2str(DFT_ratio(find(unique_elevation==q),i),2)]); %added DFTR Full -- SR 2/26/08
    end
%   print;
%   close;
end
%--------------------------------------------------------------------------
end 


%----------------------------------------------------------------------------
% Find the dftr/mfr of each individual trial
%----------------------------------------------------------------------------
trial_time = (timebin/1000):(timebin/1000):5.001;
Azims = [0 45 90 135 180 225 270 315];
Elevs = [-90 -45 0 45 90];
azimuth1 = azimuth(find(condition_num == actual_stim & fix_x == 0 & fix_y == 0));  % restricted to only 0 deg fixation trials -- CRF 3/2007
elevation1 = elevation(find(condition_num == actual_stim & fix_x == 0 & fix_y == 0));
spike_rates_dum = spike_rates(find(condition_num == actual_stim & fix_x == 0 & fix_y == 0));
temp_spike_data = data.spike_data(1,:);


% puts the binned spike count data for each trial in a cell array (all_count)
step = round(1:timebin:5002);
ii = 1;
for i = 1:5000:length(temp_spike_data)  
   tdat = temp_spike_data(i:i+4999);
   for n=1:(x_length)
       all_count{ii}(n) = sum(tdat(step(n):step(n+1)-1));
   end
   ii = ii+1;
end


% For SELECT, changed azimuth/elevation, etc. to temp_azimuth/temp_elevation to include null trials.
% This was necessary to match the trial numbers for the parsing-by-repetition in the ANOVA code below. -CRF 3/2007

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Also need to replace spike_rates_dum with the actual numbers
% from the PSTHs (2-s window shifted by delay)
temp_spike_rates2 = temp_spike_rates;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k=1:length(unique_condition_num)
    for j=1:length(unique_elevation)
        for i=1:length(unique_azimuth)          % SELECT now restricted to only 0 deg fixation trials -CRF 3/2007
            select = logical( (temp_azimuth==unique_azimuth(i)) & (temp_elevation==unique_elevation(j)) & (temp_condition_num==unique_condition_num(k)) & (temp_fix_x==0) & (temp_fix_y==0) );
            % get rid off -90 and 90 cases
            if sum(select) > 0
                act_found = find( select==1 );
                current_azi_ele = find( (Azims(plot_col) == unique_azimuth(i)) & (Elevs(plot_row) == unique_elevation(j)) );
                trial_span = span;
                for repeat=1:length(act_found)
                    
                    temp_trial_data = all_count{act_found(repeat)};
                    trial_peak = peak_time(current_azi_ele);
                    shift_trial = shift_win(current_azi_ele);

                    
                    
                    %NEW baseline test for position coding
                    pre_baseline_mean(repeat,current_azi_ele) = sum(temp_trial_data(baseline_begin:baseline_end)) / (timebin/1000);
                    post_baseline_mean(repeat,current_azi_ele) = sum(temp_trial_data(hilbert_end+5:hilbert_end+14)) / (timebin/1000);
                    baseline_diff(repeat,current_azi_ele) = post_baseline_mean(repeat,current_azi_ele) - pre_baseline_mean(repeat,current_azi_ele);

                    
                    
                    signal = temp_trial_data(baseline_begin : baseline_end);
                    temp_trial_data = temp_trial_data(1:trial_span);
                    if hilbert_end + shift_trial > trial_span
                        for ii = 1 : (hilbert_end + shift_trial - trial_span)
                            rand_sig = randperm(length(signal));
                            temp_trial_data(trial_span+ii) = signal(rand_sig(1));
                        end
                    end
                    trial_time_window = gauss_time{current_azi_ele};
                    trial_data_window = hilbert_begin+shift_trial : hilbert_end+shift_trial;

                    % corrected MFR (using delay-shifed window)
                    temp_spike_rates2(act_found(repeat)) = sum(temp_trial_data(trial_data_window)) / 2;
                    % store MFR for each dir/trial, to get stdev later
                    MFR_trial{current_azi_ele}(repeat) = sum(temp_trial_data(trial_data_window)) / 2;
                                        
                    % now we have to calculate the dft ratio for each trial
                    [trial_f trial_amp trial_phase] = FT(trial_time_window, temp_trial_data(trial_data_window), length(trial_time_window), 1, 0);
                    trial_f = round(trial_f*100)/100; % get rid of some floating point issues
                    % NOTE DFT cut-off added. CRF and SR, 3-18-2008.
                    f1_trial = mean(trial_amp(find(trial_f > 0 & trial_f <= dftcutoff)));
                    f2_trial = mean(trial_amp(find(trial_f > dftcutoff)));

                    if f2_trial == 0
                        dft_trial_all(act_found(repeat)) = 0;
                    else
                        dft_trial_all(act_found(repeat)) = f1_trial/f2_trial;
                    end

%                     % new method
%                     pol = -1;
%                     for s = find(trial_f > 0 & trial_f <= 2)
%                         amp_vel(s) = pol*trial_amp(s)*cos(trial_phase(s));
%                         amp_acc(s) = pol*trial_amp(s)*sin(trial_phase(s));
%                         pol = -pol;
%                     end
%                     v1 = mean(amp_vel(find(trial_f > 0 & trial_f <= 2)));
%                     if f2_trial == 0
%                         dft_velocity_all(act_found(repeat)) = 0;
%                     else
%                         dft_velocity_all(act_found(repeat)) = v1/f2_trial;
%                     end
%                     a1 = mean(amp_acc(find(trial_f > 0 & trial_f <= 2)));
%                     if f2_trial == 0
%                         dft_acceleration_all(act_found(repeat)) = 0;
%                     else
%                         dft_acceleration_all(act_found(repeat)) = a1/f2_trial;
%                     end

                    % phase of response
                    trial_max_c = find(trial_amp == max(trial_amp(find(trial_f > 0 & trial_f <= dftcutoff))));
                    trial_max_comp(act_found(repeat)) = trial_max_c(1);
                    trial_max_p = mean(trial_phase(trial_max_c));
                    trial_max_phase(act_found(repeat)) = trial_max_p(1);

                    if (use_max_phase == 0)
                        pol = -1;
                        for s = find(trial_f > 0 & trial_f <= dftcutoff)
                            amp_vel(s) = pol*trial_amp(s)*cos(trial_phase(s));
                            amp_acc(s) = pol*trial_amp(s)*sin(trial_phase(s));
                            if polswitch
                                pol = -pol;
                            end
                        end
                        v1 = mean(amp_vel(find(trial_f > 0 & trial_f <= dftcutoff)));
                        if f2_trial == 0
                            dft_velocity_all(act_found(repeat)) = 0;
                        else
                            dft_velocity_all(act_found(repeat)) = v1/f2_trial;
                        end
                        a1 = mean(amp_acc(find(trial_f > 0 & trial_f <= dftcutoff)));
                        if f2_trial == 0
                            dft_acceleration_all(act_found(repeat)) = 0;
                        else
                            dft_acceleration_all(act_found(repeat)) = a1/f2_trial;
                        end
                    elseif (use_max_phase == 1)
                        if polswitch
					        if trial_max_comp(act_found(repeat)) == 3
					            pol = 1;
                            else
					            pol = -1;
					        end
                        else
                            pol = -1;
                        end
                        dft_velocity_all(act_found(repeat)) = pol*cos(trial_max_phase(act_found(repeat)))*f1_trial / f2_trial;
                        dft_acceleration_all(act_found(repeat)) = pol*sin(trial_max_phase(act_found(repeat)))*f1_trial / f2_trial;
                    end
                end
                MFR_trial_mean(current_azi_ele) = mean(MFR_trial{current_azi_ele});
                MFR_trial_std(current_azi_ele) = std(MFR_trial{current_azi_ele});                
            end
        end
    end
end

% 
% 
% 
% 
% save(['Z:\Users\Chris2\tempo_backdoor\baseline_test\' files{filenum} '.mat'], 'FILE', 'baseline_diff', 'baseline_diff_avg');
% 
% 
% 
% 

% repeat for null trials only (for new significance test)
select = null_trials(select_trials) & temp_fix_x(select_trials)==0 & temp_fix_y(select_trials)==0;
if sum(select) > 0
    act_found = find( select==1 );
    current_azi_ele = 27;
    trial_span = span;
    for repeat=1:length(act_found)
        
        temp_trial_data = all_count{act_found(repeat)};
        trial_peak = peak_time(current_azi_ele);
        shift_trial = shift_win(current_azi_ele);

        signal = temp_trial_data(baseline_begin : baseline_end);
        temp_trial_data = temp_trial_data(1:trial_span);
        if hilbert_end + shift_trial > trial_span
            for ii = 1 : (hilbert_end + shift_trial - trial_span)
                rand_sig = randperm(length(signal));
                temp_trial_data(trial_span+ii) = signal(rand_sig(1));
            end
        end
        trial_time_window = gauss_time{current_azi_ele};
        trial_data_window = hilbert_begin+shift_trial : hilbert_end+shift_trial;

        % corrected MFR (using delay-shifed window)
        temp_spike_rates2(act_found(repeat)) = sum(temp_trial_data(trial_data_window)) / 2;
        % store MFR for each dir/trial, to get stdev later
        MFR_trial{current_azi_ele}(repeat) = sum(temp_trial_data(trial_data_window)) / 2;
                    
        % now we have to calculate the dft ratio for each trial
        [trial_f trial_amp trial_phase] = FT(trial_time_window, temp_trial_data(trial_data_window), length(trial_time_window), 1, 0);
        trial_f = round(trial_f*100)/100; % get rid of some floating point issues
        % NOTE DFT cut-off added. CRF and SR, 3-18-2008.
        f1_trial = mean(trial_amp(find(trial_f > 0 & trial_f <= dftcutoff)));
        f2_trial = mean(trial_amp(find(trial_f > dftcutoff)));
                    
        if f2_trial == 0
            dft_trial_all(act_found(repeat)) = 0;
        else
            dft_trial_all(act_found(repeat)) = f1_trial/f2_trial;
        end
        
%         % new method
%         pol = -1;
%         for s = find(trial_f > 0 & trial_f <= 2)
%             amp_vel(s) = pol*trial_amp(s)*cos(trial_phase(s));
%             amp_acc(s) = pol*trial_amp(s)*sin(trial_phase(s));
%             pol = -pol;
%         end
%         v1 = mean(amp_vel(find(trial_f > 0 & trial_f <= 2)));
%         if f2_trial == 0
%             dft_velocity_all(act_found(repeat)) = 0;
%         else
%             dft_velocity_all(act_found(repeat)) = v1/f2_trial;
%         end
%         a1 = mean(amp_acc(find(trial_f > 0 & trial_f <= 2)));
%         if f2_trial == 0
%             dft_acceleration_all(act_found(repeat)) = 0;
%         else
%             dft_acceleration_all(act_found(repeat)) = a1/f2_trial;
%         end
                    
        % phase of response
        trial_max_c = find(trial_amp == max(trial_amp(find(trial_f > 0 & trial_f <= dftcutoff))));
        trial_max_comp(act_found(repeat)) = trial_max_c(1);
        trial_max_p = mean(trial_phase(trial_max_c));
        trial_max_phase(act_found(repeat)) = trial_max_p(1);

        if (use_max_phase == 0)
            pol = -1;
            for s = find(trial_f > 0 & trial_f <= dftcutoff)
                amp_vel(s) = pol*trial_amp(s)*cos(trial_phase(s));
                amp_acc(s) = pol*trial_amp(s)*sin(trial_phase(s));
                if polswitch
                    pol = -pol;
                end
            end
            v1 = mean(amp_vel(find(trial_f > 0 & trial_f <= dftcutoff)));
            if f2_trial == 0
                dft_velocity_all(act_found(repeat)) = 0;
            else
                dft_velocity_all(act_found(repeat)) = v1/f2_trial;
            end
            a1 = mean(amp_acc(find(trial_f > 0 & trial_f <= dftcutoff)));
            if f2_trial == 0
                dft_acceleration_all(act_found(repeat)) = 0;
            else
                dft_acceleration_all(act_found(repeat)) = a1/f2_trial;
            end
        elseif (use_max_phase == 1)
            if polswitch
		        if trial_max_comp(act_found(repeat)) == 3
		            pol = 1;
                else
		            pol = -1;
		        end
            else
                pol = -1;
            end
            dft_velocity_all(act_found(repeat)) = pol*cos(trial_max_phase(act_found(repeat)))*f1_trial / f2_trial;
            dft_acceleration_all(act_found(repeat)) = pol*sin(trial_max_phase(act_found(repeat)))*f1_trial / f2_trial;
        end
    end
    MFR_trial_mean(current_azi_ele) = mean(MFR_trial{current_azi_ele});
    MFR_trial_std(current_azi_ele) = std(MFR_trial{current_azi_ele});   
end

% corrected MFR (using delay-shifed window)
spike_rates2 = temp_spike_rates2(~null_trials & select_trials);
spike_rates_dum2 = spike_rates2(find(condition_num == actual_stim & fix_x == 0 & fix_y == 0));


%----------------------------------------------------------------------------
% NEW SIGNIFICANCE TEST: compare single-trial DFTRs against null-trial DFTRs  --CRF 10/2007
%----------------------------------------------------------------------------
select_null = null_trials(select_trials) & temp_fix_x(select_trials)==0 & temp_fix_y(select_trials)==0;
for i = 1:26
    select = logical( (temp_azimuth==unique_azimuth(plot_col(i))) & (temp_elevation==unique_elevation(plot_row(i))) & (temp_condition_num==1) & (temp_fix_x==0) & (temp_fix_y==0) );
    if sum(select) > 0
        vector = [dft_trial_all(select) dft_trial_all(select_null)];
        group(1 : length(dft_trial_all(select))) = 1;
        group(end+1 : length(vector)) = 2;
        p_dft_vs_null(i) = anova1(vector,group,'off');
        vector = [dft_velocity_all(select) dft_velocity_all(select_null)];
        p_vel_vs_null(i) = anova1(vector,group,'off');
        vector = [dft_acceleration_all(select) dft_acceleration_all(select_null)];
        p_acc_vs_null(i) = anova1(vector,group,'off');
    end
end
p_dft_vs_null(27) = 1; p_vel_vs_null(27) = 1; p_acc_vs_null(27) = 1;

%----------------------------------------------------------------------------
% Spatial Tuning ANOVA
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
         trial_resp_matrix{k}(:,i) = [temp_spike_rates2(select_rep) NaN*ones(1,difc)];
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
dft_trial = dft_trial_all(temp_stim_type == actual_stim & temp_fix_x==0 & temp_fix_y==0 & ~null_trials & select_trials);
dft_velocity = dft_velocity_all(temp_stim_type == actual_stim & temp_fix_x==0 & temp_fix_y==0 & ~null_trials & select_trials);
dft_acceleration = dft_acceleration_all(temp_stim_type == actual_stim & temp_fix_x==0 & temp_fix_y==0 & ~null_trials & select_trials);

[DDI_mfr, var_term1] = Compute_DDI_anuk(azimuth1, elevation1, spike_rates_dum2); % added by CRF 3/2007
[DDI_full, var_term1]= Compute_DDI_anuk(azimuth1, elevation1, dft_trial); % added by CRF 3/2007
[DDI_vel, var_term1] = Compute_DDI_anuk(azimuth1, elevation1, dft_velocity);
[DDI_acc, var_term1] = Compute_DDI_anuk(azimuth1, elevation1, dft_acceleration);

% make sure they are the same length before permutation test
if length(dft_velocity) ~= length(spike_rates_dum2) | length(dft_velocity) ~= length(dft_trial)
    disp('ERROR - DFT/MFR vectors must be the same length to compute DDIs and spatial correlation');
    return
end

%permute
for i = 1:1000
    ran_ind = randperm(length(dft_velocity));
    [DDI_mfr_rand(i), var_term1] = Compute_DDI_anuk(azimuth1, elevation1, spike_rates_dum2(ran_ind));
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

%%------------------------------------------------------------------------
%% Spatial Correlation
%%------------------------------------------------------------------------
resp_mat = [];
resp_mat_vector = [];
resp_mat_std = [];
for i=1:length(unique_azimuth)
    for j=1:length(unique_elevation)
        for k=1                         % VELOCITY %
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
        for k=2                         % ACCELERATION %
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

for i=1:length(unique_azimuth)          % MEAN FIRING RATE %
    for j=1:length(unique_elevation)
        for k=3                   % changed to azimuth1/elevation1 (using spike_rates_dum2), for consistency  -CRF
            select = logical( (azimuth1==unique_azimuth(i)) & (elevation1==unique_elevation(j)));
            if (sum(select) > 0)                
                resp_mat(k, j, i) = mean(spike_rates_dum2(select));
                resp_mat_vector(k, j, i) = mean(spike_rates_dum2(select)); % for vector sum only
                for t = 1 : length(spike_rates_dum2(select));              % this is to calculate response matrix based on each trial
                    spike_temp = spike_rates_dum2(select);                 % in order to calculate error between trials in one condition
                    resp_mat_trial{k}(t, j, i) = spike_temp( t );     % t represents how many repetions each condition
                end
                resp_mat_std(k, j, i) = std(spike_rates_dum2(select));     % calculate std between trials for later DSI usage
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
% figure; contourf(resp_mat_tran_vel);
% set(gca, 'ytick', [1 2 3 4 5] );
% set(gca, 'xtick', [1 2 3 4 5 6 7 8 9] );
% set(gca, 'xticklabel','270|225|180|135|90|45|0|-45|-90');
% set(gca, 'yticklabel','90|45|0|-45|-90');
% title([FILE '  DFTR Vel']);
% colorbar;
% saveas(gcf,['C:\MATLAB6p5\work\Freq_analysis_figures\' FILE '_vel.fig']);
% close;
% 
% figure; contourf(resp_mat_tran_acc);
% set(gca, 'ytick', [1 2 3 4 5] );
% set(gca, 'xtick', [1 2 3 4 5 6 7 8 9] );
% set(gca, 'xticklabel','270|225|180|135|90|45|0|-45|-90');
% set(gca, 'yticklabel','90|45|0|-45|-90');
% title([FILE '  DFTR Acc']);
% colorbar;
% saveas(gcf,['C:\MATLAB6p5\work\Freq_analysis_figures\' FILE '_acc.fig']);
% close;
% 
% figure; contourf(resp_mat_tran_mfr);
% set(gca, 'ytick', [1 2 3 4 5] );
% set(gca, 'xtick', [1 2 3 4 5 6 7 8 9] );
% set(gca, 'xticklabel','270|225|180|135|90|45|0|-45|-90');
% set(gca, 'yticklabel','90|45|0|-45|-90');
% title([FILE '  MFR']);
% colorbar;
% saveas(gcf,['C:\MATLAB6p5\work\Freq_analysis_figures\' FILE '_mfr.fig']);
% close;
% % ----------------------------------------------------------------------------

%Correlate the mean firing rate contour plot with the dftr contour plot.

% ---------------
% THIS IS USING AVERAGED SINGLE-TRIAL DFTR COMPONENTS (too noisy?)
% ---------------
% count = 1;
% for i = 1:length(unique_azimuth)
%     for j = 1:length(unique_elevation)
%         if (unique_elevation(j) == -90 | unique_elevation(j) == 90) & unique_azimuth(i)~= 0 
%             k = 1;
%         else
%             vel_data(count) = resp_mat_tran(1,j,i);
%             acc_data(count) = resp_mat_tran(2,j,i);
%             mfr_data(count) = resp_mat_tran(3,j,i);
%             count = count+1;
%         end
%     end
% end

% %---------------
% % THIS IS USING DFTR COMPONENTS COMPUTED FROM AVERAGED PSTHs (better)
% %---------------
count = 1;
for i = 1:length(unique_azimuth)
    for j = 1:length(unique_elevation)
        if (unique_elevation(j) == -90 | unique_elevation(j) == 90) & unique_azimuth(i)~= 0
            k = 1;  % that's weird..
        else
            vel_data(count) = vel_sum(j,i);
            acc_data(count) = acc_sum(j,i);
            mfr_data(count) = resp2{current_stim}(j,i);
            count = count + 1;
        end
    end
end

% % Flags the nonsig cells with 999 rather than a (probably) meaningless correlation/P value
% vel_mfr = 999;
% p_vel_mfr = 999;
% acc_mfr = 999;
% p_acc_mfr = 999;
% acc_vel = 999;
% p_acc_vel = 999;
%
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
% foldername = ('C:\MATLAB6p5\work\Frequency_Analysis\');
% foldername = ('C:\Aihua\z_TempOutputs\Frequency_Analysis\');
foldername = ('Z:\Users\Tanya\Frequency_Analysis_output\');

outfile1 = [foldername 'all_data_transvest_Priam1.dat'];
outfile2 = [foldername 'all_data_transvest_Priam2.dat'];

for i = 1:27
    sprint_txt = ['%s\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t'];
    if i == 27
        buff = sprintf(sprint_txt, FILE, current_stim, i, 99, 99, fourier_ratio(i), dft_p(i), fourier_ratio_vel(i), dft_p_vel(i), fourier_ratio_acc(i), dft_p_acc(i), max_phase(i), shift_win(i)*(timebin/1000), p_dft_vs_null(i), p_vel_vs_null(i), p_acc_vs_null(i), bpr(i), vel_pct(i), acc_pct(i), vel_pct_adj(i), acc_pct_adj(i), emptyerror1(i), emptyerror2(i), MFR_trial_mean(i), MFR_trial_std(i), repetitions, dft_p_vs_Null(i));
    else
        buff = sprintf(sprint_txt, FILE, current_stim, i, unique_azimuth(plot_col(i)), unique_elevation(plot_row(i)), fourier_ratio(i), dft_p(i), fourier_ratio_vel(i), dft_p_vel(i), fourier_ratio_acc(i), dft_p_acc(i), max_phase(i), shift_win(i)*(timebin/1000), p_dft_vs_null(i), p_vel_vs_null(i), p_acc_vs_null(i), bpr(i), vel_pct(i), acc_pct(i), vel_pct_adj(i), acc_pct_adj(i), emptyerror1(i), emptyerror2(i), MFR_trial_mean(i), MFR_trial_std(i), repetitions, dft_p_vs_Null(i));
    end
    fid = fopen(outfile1, 'a');
    fprintf(fid, '%s', buff);
    fprintf(fid, '\r\n');
    fclose(fid);
end

sprint_txt = ['%s\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t'];
buff = sprintf(sprint_txt, FILE, current_stim, P_anova(current_stim), P_anova_full(current_stim), P_anova_vel(current_stim), P_anova_acc(current_stim), DDI_mfr, DDI_full, DDI_vel, DDI_acc, p_ddi_mfr, p_ddi_full, p_ddi_vel, p_ddi_acc, vel_mfr, p_vel_mfr,...
    acc_mfr, p_acc_mfr, acc_vel, p_acc_vel, mfr_pref_az, mfr_pref_el, mfr_pref_amp, pref_az_vel, pref_el_vel, pref_amp_vel, pref_az_acc, pref_el_acc, pref_amp_acc);  
fid = fopen(outfile2, 'a');
fprintf(fid, '%s', buff);
fprintf(fid, '\r\n');
fclose(fid);
% 	for i = 1:27
%          sprint_txt = ['%s\t %f\t %f\t %f\t'];
%          if i == 27
%              buff = sprintf(sprint_txt, FILE, i, bpr(i));  
%          else
%              buff = sprintf(sprint_txt, FILE, i, bpr(i));  
%          end
%          fid = fopen(outfile1, 'a');
%          fprintf(fid, '%s', buff);
%          fprintf(fid, '\r\n');
%          fclose(fid);
%      end

time = clock; disp(['End Time = ' num2str([time(4) time(5) time(6)])]);

return;