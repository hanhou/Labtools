%----------------------------------------------------------------------------------------------------------------------
%Frequency_Analysis.m calculates the DFTR, phase, latency, velocity and acceleration
%components, preferred directions (for MFR, velocity and acceleration),
%DDI's, vector sum amplitudes, correlation of 3d tuning and anova for neuronal responses

function Frequency_Analysis_alternate(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);

time = clock; disp(['Start Time = ' num2str([time(4) time(5) time(6)])]);

Path_Defs;
ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP

%get the column of values for azimuth and elevation and stim_type
temp_azimuth = data.moog_params(AZIMUTH,:,MOOG);
temp_elevation = data.moog_params(ELEVATION,:,MOOG);
temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG);
temp_amplitude = data.moog_params(AMPLITUDE,:,MOOG);

%now, get the firing rates for all the trials
temp_spike_rates = data.spike_rates(SpikeChan, :);                                                                                                                             

%get indices of any NULL conditions (for measuring spontaneous activity
null_trials = logical( (temp_azimuth == data.one_time_params(NULL_VALUE)) );

%now, remove trials from direction and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:length(temp_azimuth);		% a vector of trial indices
bad_tri = find(temp_spike_rates > 3000);   % cut off 3k frequency which definately is not cell's firing response
if ( bad_tri ~= NaN)
   select_trials= ( (trials >= BegTrial) & (trials <= EndTrial) & (trials~=bad_tri) );
else 
   select_trials= ( (trials >= BegTrial) & (trials <= EndTrial) ); 
end

azimuth = temp_azimuth(~null_trials & select_trials);
elevation = temp_elevation(~null_trials & select_trials);
stim_type = temp_stim_type(~null_trials & select_trials);
amplitude = temp_amplitude(~null_trials & select_trials);
spike_rates = temp_spike_rates(~null_trials & select_trials);

unique_azimuth = munique(azimuth');
unique_elevation = munique(elevation');
unique_stim_type = munique(stim_type');
unique_amplitude = munique(amplitude');

resp_mat = [];
resp_mat_std = [];
for i=1:length(unique_azimuth)
    for k=1: length(unique_stim_type)
        select = logical( (azimuth==unique_azimuth(i)) & (elevation==0) & (stim_type==unique_stim_type(k)) );
        if sum(select) > 0
            resp_mat(i,k) = mean(spike_rates(select));
            resp_mat_std(i,k) = std(spike_rates(select));
        else
            resp_mat(i,k) = NaN;
            resp_mat_std(i,k) = NaN;
        end
    end
end

trials_per_rep = length(unique_azimuth)*length(unique_elevation)*length(unique_stim_type)+1;
num_reps = floor( (EndTrial-(BegTrial-1)) / trials_per_rep); 

save(['Z:\Users\Chris2\tempo_backdoor\' FILE '.mat'], 'Protocol', 'unique_stim_type', 'unique_azimuth', 'resp_mat', 'resp_mat_std', 'num_reps');

% FISHER = slope^2 / variance
% 
% function Frequency_Analysis_alternate(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
% 
% disp(' '); time = clock; disp(['Start Time = ' num2str([time(4) time(5) time(6)])]);
% 
% Path_Defs;
% ProtocolDefs;  %contains protocol specific keywords - 1/4/01 BJP
% 
% plot_psth = 0;
% plot_spatial = 0;
% 
% dftcutoff = 2.5;    % Select the DFTR cut-off 
% polswitch = 1;      % Select 1 for switching signs on the component phases
% use_max_phase = 0;
% 
% % load('Z:\LabTools\Matlab\TEMPO_Analysis\ProtocolSpecific\MOOG\time_vel_acc.mat');
% % load('Z:\Data\Tempo\Batch Files\Suhrud\time_vel_acc.mat');
% % load('C:\MATLAB6p5\work\time_vel_acc.mat');
% load('Z:\Users\Chris2\time_vel_acc.mat');
% 
% % load('Z:\LabTools\Matlab\TEMPO_Analysis\ProtocolSpecific\MOOG\bpr_lookup_table.mat');
% % load('Z:\Data\Tempo\Batch Files\Suhrud\bpr_lookup_table.mat');
% 
% if Protocol == 100 || Protocol == 104
% 	%get the column of values for azimuth and elevation and stim_type
% 	temp_azimuth = data.moog_params(AZIMUTH,:,MOOG);
% 	temp_elevation = data.moog_params(ELEVATION,:,MOOG);
% 	temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG);
% 	temp_spike_data = data.spike_data(SpikeChan,:);
% 	% fixation data added by CRF, 3/2007
% 	temp_fix_x    =  data.moog_params(FIX_X,:,MOOG);
% 	temp_fix_y    =  data.moog_params(FIX_Y,:,MOOG);
% 	temp_fix_x(isnan(temp_fix_x)) = 0;
% 	temp_fix_y(isnan(temp_fix_y)) = 0;
% elseif Protocol == 112
% 	%get the column of values for azimuth and elevation and stim_type
% 	temp_azimuth = data.moog_params(ROT_AZIMUTH,:,MOOG);
% 	temp_elevation = data.moog_params(ROT_ELEVATION,:,MOOG);
% 	temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG);
% 	temp_spike_data = data.spike_data(SpikeChan,:);
% 	% fixation data added by CRF, 3/2007
% 	temp_fix_x    =  data.moog_params(FIX_X,:,MOOG);
% 	temp_fix_y    =  data.moog_params(FIX_Y,:,MOOG);
% 	temp_fix_x(isnan(temp_fix_x)) = 0;
% 	temp_fix_y(isnan(temp_fix_y)) = 0;    
% end
%     
%     
% %mean firing rate of each trial depending on the start and stop offsets
% temp_spike_rates = data.spike_rates(SpikeChan, :);
% 
% %--------------------------------------------------------------
% % Some cells have more than the usual 26 directions, so we need to remove them
% extra_directions =[];
% if length(munique(temp_azimuth')) > 9  % it's 9 because the null trials (-9999) are included
%     unique_azimuth_withnull = [0;45;90;135;180;225;270;315;-9999];
%     unique_elevation_withnull = [-90;-45;0;45;90;-9999];
%     extra_azis = ones(1,length(temp_azimuth));
%     extra_eles = ones(1,length(temp_elevation));
%     for n = 1:length(unique_azimuth_withnull)
%         extra_azis(find(temp_azimuth == unique_azimuth_withnull(n))) = 0;
%     end
%     for m = 1:length(unique_elevation_withnull)
%         extra_eles(find(temp_elevation == unique_elevation_withnull(m))) = 0;
%     end
%     extra_directions = (extra_azis | extra_eles);
%     temp_elevation(extra_directions) = [];
%     temp_stim_type(extra_directions) = [];
%     temp_fix_x(extra_directions) = [];
%     temp_fix_y(extra_directions) = [];
%     temp_spike_rates(extra_directions) = [];    
%     temp_azimuth(extra_directions) = [];
%     
%     EndTrial = EndTrial - sum(extra_directions(1:EndTrial-(BegTrial-1)));
% end
% %--------------------------------------------------------------
% 
% %get indices of any NULL conditions (for measuring spontaneous activity)
% null_trials = logical( (temp_azimuth == data.one_time_params(NULL_VALUE)) );
% 
% %now, remove trials from direction and spike_rates that do not fall between BegTrial and EndTrial
% trials = 1:length(temp_azimuth); % a vector of trial indices
% %at any given point, cell response cannot be greater than 1
% abnormal = find(temp_spike_data > 1);
% temp_spike_data(1,abnormal) = 1;
% bad_trials = find(temp_spike_rates > 3000);   % cut off 3k frequency which definately is not cell's firing response
% 
% if ( bad_trials ~= NaN)
%     select_trials= ( (trials >= BegTrial) & (trials <= EndTrial) & (trials~=bad_trials) );
% else
%     select_trials= ( (trials >= BegTrial) & (trials <= EndTrial) );
% end
% 
% azimuth = temp_azimuth(~null_trials & select_trials);
% elevation = temp_elevation(~null_trials & select_trials);
% stim_type = temp_stim_type(~null_trials & select_trials);
% fix_x = temp_fix_x(~null_trials & select_trials);
% fix_y = temp_fix_y(~null_trials & select_trials);
% spike_rates = temp_spike_rates(~null_trials & select_trials);
% 
% fix_x_withnull = temp_fix_x(select_trials);
% fix_y_withnull = temp_fix_y(select_trials);
% spike_rates_withnull = temp_spike_rates(select_trials);
% 
% unique_azimuth = munique(azimuth');
% unique_elevation = munique(elevation');
% unique_stim_type = munique(stim_type');
% unique_fix_x    =  munique(fix_x');
% unique_fix_y    =  munique(fix_y');
% 
% temp_condition_num = temp_stim_type;
% condition_num = stim_type;
% h_title{1}='Vestibular';
% h_title{2}='Visual';
% h_title{3}='Combined';
% unique_condition_num = munique(condition_num');
% 
% % --------------------------------------------
% % select bin size -- 50 ms (20 Hz) is default
% 
% % timebin = 50;
% % time = time50; vel = vel50; acc = acc50;
% 
% timebin = 20;
% time = time20; vel = vel20; acc = acc20;
% 
% % timebin = 16.6667;
% % time = time16; vel = vel16; acc = acc16;
% 
% % timebin = 1;
% % time = time1; vel = vel1; acc = acc1;
% % --------------------------------------------
% 
% % sample frequency depends on test duration
% % frequency=length(temp_spike_data)/length(select_trials); % actually its 5000 now for 5 sec 
% frequency = 5000;
% 
% % length of x-axis
% x_length = round(frequency/timebin);
% % x-axis for plot PSTH
% x_time=1:x_length;
% x_time_bincenter = timebin/2000;            % using bin-centers now: start at time = 0, plus a half-bin...
% for nn = 1 : round(5/(timebin/1000)) - 1    % then fill it out with a full 2-seconds (i.e., 40 bins with timebin = 50 ms)
%     x_time_bincenter(end+1) = x_time_bincenter(end) + (timebin/1000);
% end
% 
% % remove 'extra direction' trials from spikedata
% discard_trials = find(extra_directions);
% for i = 1 : length(discard_trials)
%     temp_spike_data( 1, ((discard_trials(i)-1)*frequency +1) : discard_trials(i)*frequency ) = 99;
% end
% temp_spike_data = temp_spike_data( temp_spike_data~=99 );
% 
% % remove trials outside Begtrial~Endtrial
% discard_trials = [find(trials <BegTrial | trials >EndTrial)];
% for i = 1 : length(discard_trials)
%     temp_spike_data( 1, ((discard_trials(i)-1)*frequency +1) : discard_trials(i)*frequency ) = 99;
% end
% spike_data_withnull = temp_spike_data( temp_spike_data~=99 );
% 
% % remove null trials and trials outside Begtrial~Endtrial
% discard_trials = [find(null_trials==1 | trials <BegTrial | trials >EndTrial)];
% for i = 1 : length(discard_trials)
%     temp_spike_data( 1, ((discard_trials(i)-1)*frequency +1) :  discard_trials(i)*frequency ) = 99;
% end
% spike_data = temp_spike_data( temp_spike_data~=99 );
% 
% % count spikes from raster data (spike_data)
% max_count = 1;
% time_step = 1;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% current_stim = 1;  % vestibular=1, visual=2
% actual_stim = current_stim;
% 
% % For visual-only blocks, the index 'k' (below) will only reach 1, and no data
% % will be created at index k = 2.  Therefore must change current_stim to 1. -CRF
% if unique_condition_num == 2 & current_stim == 2
%     current_stim = 1;
%     actual_stim = 2;
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% for k=1:length(unique_condition_num)  % it's the stim_type
%     for j=1:length(unique_elevation)
%         for i=1:length(unique_azimuth)         % SELECT now restricted to only 0 deg fixation trials -- CRF 3/2007
%             select = logical(azimuth==unique_azimuth(i) & elevation==unique_elevation(j) & condition_num==unique_condition_num(k) & fix_x==0 & fix_y==0 );            
%             if sum(select) > 0
%                 resp{k}(j,i) = mean(spike_rates(select));
%                 act_found = find( select==1 );
%                 % count spikes per timebin on every same condition trials
%                 clear temp_count dummy_count;
%                 for repeat=1:length(act_found) 
%                     for n=1:(x_length)
%                         temp_count(repeat,n)=sum( spike_data( 1, frequency*(act_found(repeat)-1) + time_step : round(frequency*(act_found(repeat)-1)+n*timebin) ) );
%                         dummy_count{repeat}(n) = sum(spike_data(1,(frequency*(act_found(repeat)-1)+time_step):round(frequency*(act_found(repeat)-1)+n*timebin)));
%                         if timebin == 16.6667
%                             if floor((n+1)/3) == (n+1)/3
%                                 time_step = floor(time_step + timebin);
% 						    else                               % kluge for 16.6667 ms timebin (just rounding was missing bins in spike_data)
%                                 time_step = round(time_step + timebin);
% 							end
%                         else
%                             time_step = time_step + timebin;   % this may need to be modified for timebins other than 50
%                         end
%                     end
%                     time_step=1;
%                     
%                     if k == current_stim
%                         count_condition{i,j,repeat} = dummy_count{repeat};
%                     end
%                 end
%                 
%                 count_y_trial{i,j,k}(:,:) = temp_count;  % each trial's PSTH in vestibular condition
%                 if k == 1
%                     count_y_trial1{i,j}(:,:) = temp_count;
%                 end
%                 
%                 dim=size(temp_count);
%                 if dim(1) > 1;
%                     count_y{i,j,k} = mean(temp_count);
%                 else
%                     count_y{i,j,k}= temp_count;     % for only one repetition cases
%                 end
%             else
%                 resp{k}(j,i) = 0; 
%                 count_y{i,j,k}=0;
%             end   
%             % normalize count_y
%             if max(count_y{i,j,k})~=0;
%                 count_y_norm{i,j,k}=count_y{i,j,k} / max(count_y{i,j,k});
%             else
%                 count_y_norm{i,j,k}=0;
%             end
%         end
%     end
%     % now find the peak
%     [row_max, col_max] = find( resp{k}(:,:)==max(max(resp{k}(:,:))) );
%     % it is likely there are two peaks with same magnitude, choose the first one arbitraly
%     row_m{k}=row_max(1);
%     col_m{k}=col_max(1);
%     if max(count_y{col_max(1), row_max(1), k})~=0;
%         count_y_max{k} = count_y{col_max(1), row_max(1), k} / max(count_y{col_max(1), row_max(1), k});
%     else
%         count_y_max{k} =0;
%     end
%     % find the largest y to set scale later
%     if max(count_y{col_max(1), row_max(1), k}) > max_count
%         max_count = max(count_y{col_max(1), row_max(1), k});
%     end    
% end
% 
% % Do the same, but for null trials (for new significance test on DFTR below) -- CRF 10/2007
% select = null_trials(select_trials) & temp_fix_x(select_trials)==0 & temp_fix_y(select_trials)==0;
% if sum(select) > 0
%     resp_null = mean(spike_rates_withnull(select));
%     act_found = find( select==1 );
%     % count spikes per timebin on every same condition trials
%     clear temp_count dummy_count;
%     for repeat=1:length(act_found) 
%         for n=1:(x_length)
%             temp_count(repeat,n)=sum( spike_data_withnull( 1, frequency*(act_found(repeat)-1) + time_step : round(frequency*(act_found(repeat)-1)+n*timebin) ) );
%             dummy_count{repeat}(n) = sum(spike_data_withnull(1,(frequency*(act_found(repeat)-1)+time_step):round(frequency*(act_found(repeat)-1)+n*timebin)));
%             if timebin == 16.6667
%                 if floor((n+1)/3) == (n+1)/3
%                     time_step = floor(time_step + timebin);
% 			    else                               % kluge for 16.6667 ms timebin (just rounding was missing bins in spike_data)
%                     time_step = round(time_step + timebin);
% 				end
%             else
%                 time_step = time_step + timebin;   % this may need to be modified for timebins other than 50
%             end
%         end
%         time_step=1;
%     end
%     count_y_trial_null = temp_count;
%     dim=size(temp_count);
%     if dim(1) > 1;
%         count_y_null = mean(temp_count);
%     else
%         count_y_null = temp_count;     % for only one repetition cases
%     end
% else
%     resp_null=0; 
%     count_y_null=0;
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % order in which the directions are plotted
% plot_col = [1 1 1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 6 6 6 7 7 7 8 8 8];
% plot_row = [5 4 3 2 1 4 3 2 4 3 2 4 3 2 4 3 2 4 3 2 4 3 2 4 3 2];
% 
% % find when in the 5-sec trial period that the spike data cuts off
% sum_dat = zeros(1,x_length);
% for j = 1:26    
%     sum_dat = sum_dat + count_y{plot_col(j), plot_row(j), current_stim};
% end
% % and accordingly set length of data to read from PSTH
% for ii = 1:x_length
%     if sum(sum_dat(ii:end)) == 0
%         span = ii;
%         break
%     end
% end
% 
% if span < .80*x_length
%     span = round(.68*x_length); % usually this is how much data there is
% else % this is one of katsu's cells
%     span = round(.80*x_length);
% end
% 
% % ----------------------------------------------------------------------------------------------
% % Define some time points for segmenting the response into the segment of interest -CRF 6/2007
% % ----------------------------------------------------------------------------------------------
% 
% % The stimulus velocity peaks at time t = 2092 ms, however the 'onset' is actually t = 996 ms, according to the event code* (4).
% % Nevertheless, for simplicity we will assume a 2-second stimulus window beginning at t = 1092 ms and ending at t = 3092.
% % *(Oddly, the stim offset code (5) occurs at t = 3006, so there's an extra 10 ms in there somehow.)
% 
% % Take baseline data as the 200 ms before stimulus onset (t = 996 ms)
% tempbins = find(x_time_bincenter > .796);
% baseline_begin = tempbins(1);
% tempbins = find(x_time_bincenter < .996);
% baseline_end = tempbins(end);
% 
% % Data for hilbert transform (temporal envelope) can be the full 2-second window
% tempbins = find(x_time_bincenter > 1.092);
% hilbert_begin = tempbins(1);
% tempbins = find(x_time_bincenter < 3.092);
% hilbert_end = tempbins(end);
% 
% % 
% % 
% % % TEMP: expanded window to address 'edge effect' concern
% % hilbert_begin = hilbert_begin - (span-hilbert_end);
% % hilbert_end = span;
% % 
% % 
% 
% % Data from 1400 to 2700 ms into stim period, from which the peak is taken
% % (but only as a failsafe for when center of mass doesn't fall in this range, which is relatively rare)
% tempbins = find(x_time_bincenter > 1.492);
% peak_begin = tempbins(1);
% tempbins = find(x_time_bincenter < 2.792);
% peak_end = tempbins(end);
% 
% % Finally, find the bin closest to the stimulus peak velocity (2.092 s)
% peak_stim = find(abs(x_time_bincenter - 2.092) == min(abs(x_time_bincenter - 2.092)));
% peak_stim = peak_stim(1);
% 
% for j = 1:27
%     emptyerror1(j) = 0;
%     emptyerror2(j) = 0;
% end
% 
% 
% % -----------------------------------------------------------------------------------------------------------
% % -----------------------------------------------------------------------------------------------------------
% % -----------------------------------------------------------------------------------------------------------
% % -----------------------------------------------------------------------------------------------------------
% % perform frequency analysis---------------------------------------------------------------------------------
% % -----------------------------------------------------------------------------------------------------------
% % -----------------------------------------------------------------------------------------------------------
% % -----------------------------------------------------------------------------------------------------------
% % -----------------------------------------------------------------------------------------------------------
% for j=1:27
%     if j == 27
%         gdat = count_y_null;
%     else
%         gdat = count_y{plot_col(j), plot_row(j), current_stim};
%     end
%    
%     %     NEW: baseline diff post-pre (see single trials loop below) to plot spatial tuning of position signal
%     pre_baseline_mean_avg(j) = sum(gdat(baseline_begin:baseline_end)) / (timebin/1000);
%     post_baseline_mean_avg(j) = sum(gdat(hilbert_end+5:hilbert_end+14)) / (timebin/1000);  % hard coded for now
%     baseline_diff_avg(j) = post_baseline_mean_avg(j) - pre_baseline_mean_avg(j);
%         
%     % %  %%NEW STUFF APR 15th BEGIN  
%     datatogetmean = gdat;
%     meanfiring = mean(datatogetmean(30:60));
%     
%     sizeofbins = 10;
%     slider = 1;
%     startat = hilbert_begin;
%     nn = 1;
%     endat = startat +sizeofbins-1;
%     while endat < hilbert_end
%         if sum(datatogetmean(startat : endat) )~=0 
%             meandata(nn) = mean(datatogetmean(startat : endat) );
%         else 
%             meandata(nn) = 0;
%         end
%         startat = startat + slider;
%         endat = startat + sizeofbins -1;
%         nn = nn + 1;
%     end
%     
%     maxmean = meandata(find(meandata == max(meandata)));
%     bpr(j) = meanfiring/maxmean(1);
%     
% % % %     figure(j+1); bar(x_time, gdat);
% % % %     hold on; plot(x_time(slider*find(meandata == max(meandata)))+hilbert_begin, maxmean, 'r*');
% % % end
% % %     
% % % %%NEW STUFF APR 15 END
% % %     
%     total_data{j} = gdat(1:span);
%     %---------------------------------------
%     % find hilbert transform to get peak.
%     shift_win(j) = 0;
%     flagged_dir(j) = 0;
%     % running mean
%     binsize = 3; % take mean of 3 neighbour bins
%     slidesize = 1; % with a slide of 1 bin
%     nn = 1;
%     start_bin = 1;
%     end_bin = start_bin + binsize -1;
%     while end_bin < x_length % the length of gdat
%         if sum( gdat(start_bin : end_bin) )~=0 % avoid divided by 0 later
%             runmean(nn) = mean( gdat(start_bin : end_bin) );
%         else 
%             runmean(nn) = 0;
%         end
%         start_bin = start_bin + slidesize;
%         end_bin = start_bin + binsize -1;
%         nn = nn + 1;
%     end
%     hildata = zeros(1,x_length);
%     hildata(2:length(runmean)+1) = runmean;
%     dummy_dat{j} = hildata(hilbert_begin : hilbert_end) - mean(hildata(baseline_begin : baseline_end));
%     x_time_bincenter_hil = x_time_bincenter(hilbert_begin : hilbert_end);
%     hil{j} = hilbert(dummy_dat{j});
%     abb = abs(hil{j});
% %     HT = abb - mean(abb(round(.01*x_length) : round(.03*x_length))); % remove the mean before finding the center of mass or else it will be biased (???)
%     HT = abb - min(abb);
%     
%     % Find four largest bins as candidates for the peak of the envelope.
%     % The bin with the largest mean over 3 bins will be the peak.
%     % This is only used when the center of mass falls outside the range 1.5s to 2.7s (not often)
%     dum_HT = HT(peak_begin-hilbert_begin+1 : peak_end-hilbert_begin+1);
%     mm = find(dum_HT == max(dum_HT));
%     m(1) = mm(1);
%     dum_HT(m(1)) = 0;
%     mm = find(dum_HT == max(dum_HT));
%     m(2) = mm(1);
%     dum_HT(m(2)) = 0;
%     mm = find(dum_HT == max(dum_HT));
%     m(3) = mm(1);
%     dum_HT(m(3)) = 0;
%     mm = find(dum_HT == max(dum_HT));
%     m(4) = mm(1);
%     dum_HT(m(4)) = 0;
%     dum_HT = HT(peak_begin-hilbert_begin+1 : peak_end-hilbert_begin+1);
%     for i = 1:4
%         if m(i) == length(dum_HT)
%             me(i) = mean(dum_HT(m(i)-1:m(i)));
%         elseif m(i) == 1;
%             me(i) = mean(dum_HT(m(i):m(i)+1));
%         else
%             me(i) = mean(dum_HT(m(i)-1:m(i)+1));
%         end
%     end
%     maxbin = find(me == max(me));
%     peak_temp = m(maxbin) + peak_begin - 1;
% 
%     H_time = x_time_bincenter(hilbert_begin : hilbert_end);
%     CMsum = HT .* H_time;
%     if sum(HT) > 0
%         center_mass = sum(CMsum)/sum(HT);
%         peak_t = find( abs(x_time_bincenter - center_mass) == min(abs(x_time_bincenter - center_mass)) );% point which is closest to the center of mass is used.
%         if x_time_bincenter(peak_t) >= 1.4 & x_time_bincenter(peak_t) <= 2.7
%             peak_time(j) = peak_t(1);
%         else 
%             peak_time(j) = peak_temp(1);
%             flagged_dir(j) = 1;
%         end
%     else
%         peak_time(j) = peak_stim;
%     end
% 
%     shift_win(j) = peak_time(j) - peak_stim; % 'delay' (shift of window for fourier transform to get rid of phase shift caused by delayed response)
% %     shift_win(j) = 10;  % temporarily fix delay
%     
%     %if the window is being shifted to extend into the region where the
%     %recording has stopped, then we will randomly sample from the front of
%     %the response and fill out the window.
%     signal = total_data{j}(baseline_begin : baseline_end);
%     span_new(j) = span;
%     
%     if hilbert_end + shift_win(j) > span
%         for i = 1:(hilbert_end + shift_win(j) - span)
%             rand_sig = randperm(length(signal));
%             total_data{j}(span+i) = signal(rand_sig(1));
%             span_new(j) = span+i;
%         end
%     end
%     
%     % window is shifted depending on the delay of the peak of the envelope (from the hilbert trans), relative to the peak of the stimulus
%     gauss_dat{j} = total_data{j}(hilbert_begin+shift_win(j) : hilbert_end+shift_win(j));
% %     gauss_dat{j} = total_data{j}(hilbert_begin+shift_win(j)-1 : hilbert_end+shift_win(j)-1);      %******************************************
%     gauss_time{j} = x_time_bincenter(hilbert_begin+shift_win(j) : hilbert_end+shift_win(j));        % NOT SURE HERE, BECAUSE SEGMENT HAS 40 BINS, AND
% %     gauss_time{j} = x_time_bincenter(hilbert_begin+shift_win(j)-1 : hilbert_end+shift_win(j)-1);  % THUS CANNOT BE CENTERED ON ONE PARTICULAR BIN.
%                                                                                                     % FOR 50-MS TIMEBIN, THIS CAUSES A PHASE DIFF OF 9 DEG.
%                                                                                                     % ******************************************
%     % Also recompute MFR with the delay-shifted window
%     if j == 27
%         resp2_null = sum(total_data{j}(hilbert_begin+shift_win(j) : hilbert_end+shift_win(j))) / 2;
%     else
%         resp2{current_stim}(plot_row(j),plot_col(j)) = sum(total_data{j}(hilbert_begin+shift_win(j) : hilbert_end+shift_win(j))) / 2;
%     end
%                                                                                                     
% %     time_short = time(hilbert_begin:hilbert_end);
% %     vel_short = vel(hilbert_begin:hilbert_end);
% %     acc_short = acc(hilbert_begin:hilbert_end);
% 
%     %calculate DFT ratio
%     %40 pt DFT with mean removed
%     [f, amp, resp_phase] = FT(gauss_time{j}, gauss_dat{j}, length(gauss_time{j}), 1, 0);
%     f = round(f*100)/100; % get rid of some floating point issues
%     
%     % NOTE DFT cut-off added. CRF and SR, 3-18-2008.
%     
%     f1 = mean(amp(find(f > 0 & f <= dftcutoff)));
%     f2 = mean(amp(find(f > dftcutoff)));
%     if f2 == 0
%         fourier_ratio(j) = 0;
%     else        %DFTR
%         fourier_ratio(j) = f1/f2;
%     end
% 
%     %phase of response
%     max_c = find(amp == max(amp(find(f > 0 & f <= dftcutoff))));
%     max_comp(j) = max_c(1);
%     max_p = mean(resp_phase(max_c));
%     max_phase(j) = max_p(1);
% 
%     %  New approach to vel and accel DFTR: multiply the amplitude of each frequency component by the 
% 	%  cosine (vel) or sine (acc) of its phase angle (essentially scaling each by a measure of how
% 	%  close its phase is to the phase of the vel and accel stimulus profiles), then recompute a separate 
% 	%  vel-DFTR and acc-DFTR from these scaled amplitudes.     GCD and CRF, 7-20-2007
% 
%     if (use_max_phase == 0)
%         pol = -1;
%         for s = find(f > 0 & f <= dftcutoff)
%             amp_vel(s) = pol*amp(s)*cos(resp_phase(s));
%             amp_acc(s) = pol*amp(s)*sin(resp_phase(s));
%             if polswitch
%                 pol = -pol;
%             end
%         end
%         v1 = mean(amp_vel(find(f > 0 & f <= dftcutoff)));
%         if f2 == 0
%             fourier_ratio_vel(j) = 0;
%         else
%             fourier_ratio_vel(j) = v1/f2;
%         end
%         a1 = mean(amp_acc(find(f > 0 & f <= dftcutoff)));
%         if f2 == 0
%             fourier_ratio_acc(j) = 0;
%         else
%             fourier_ratio_acc(j) = a1/f2;
%         end
%     elseif (use_max_phase == 1)
%         if polswitch
% 	        if max_comp(j) == 3
% 	            pol = 1;
%             else
% 	            pol = -1;
% 	        end
%         else
%             pol = -1;
%         end
%         fourier_ratio_vel(j) = pol*cos(max_phase(j))*f1 / f2;
%         fourier_ratio_acc(j) = pol*sin(max_phase(j))*f1 / f2;
%     end
%     
%     % Define a vel and acc 'percentage'
%     dft_sum(j) = abs(fourier_ratio_vel(j)) + abs(fourier_ratio_acc(j));
%     vel_pct(j) = fourier_ratio_vel(j) / dft_sum(j);
%     vel_pct(j) = round(vel_pct(j)*100) / 100;
%     acc_pct(j) = fourier_ratio_acc(j) / dft_sum(j);
%     acc_pct(j) = round(acc_pct(j)*100) / 100;
% 
% %     % Adjust vel/acc percentages using a lookup-table approach, based on fake data
% %     bpr(j) = round(bpr(j)*100) / 100;
% %     vel_pct_adj_set = [];
% %     acc_pct_adj_set = [];
% %     emptyerror1(j) = 0;
% %     emptyerror2(j) = 0;
% %     
% %     if isempty(find(bpr_lookup == bpr(j) & vel_pct_output == vel_pct(j)))
% %         emptyerror1(j) = 1;
% %         vel_pct_adj_set = vel_pct(j);
% %     else
% %         vel_pct_adj_set = vel_pct_input( find(bpr_lookup == bpr(j) & vel_pct_output == vel_pct(j)) );
% % 	end
% % 	
% %     if isempty(find(bpr_lookup == bpr(j) & acc_pct_output == acc_pct(j)))
% %         emptyerror2(j) = 1;
% %         acc_pct_adj_set = acc_pct(j);
% %     else
% %         acc_pct_adj_set = acc_pct_input( find(bpr_lookup == bpr(j) & acc_pct_output == acc_pct(j)) );
% % 	end
% %     
% %     index = 1;
% %     for a = 1:length(vel_pct_adj_set)
% %         for b = 1:length(acc_pct_adj_set)
% %             if abs(vel_pct_adj_set(a)) + abs(acc_pct_adj_set(b)) == 1
% %                 vel_pct_good(index) = vel_pct_adj_set(a);
% %                 acc_pct_good(index) = acc_pct_adj_set(b);
% %                 index = index + 1;
% %             end
% %         end
% %     end
% % 
% %     if emptyerror1(j)
% %         vel_pct_adj(j) = vel_pct(j);
% %     else
% %         vel_pct_adj(j) = mean(unique(vel_pct_good));
% %     end
% %     if emptyerror2(j)
% %         acc_pct_adj(j) = acc_pct(j);
% %     else
% %         acc_pct_adj(j) = mean(unique(acc_pct_good));
% %     end
% %     
% %     % Lastly, recompute vel and acc DFTRs via the new _pcts
% %     vel_dftr_adj(j) = vel_pct_adj(j) * dft_sum(j);
% %     acc_dftr_adj(j) = acc_pct_adj(j) * dft_sum(j);
% 
% 
% % % % don't adjust dftr's for now
% % %     
% % %     fourier_ratio_vel(j) = vel_dftr_adj(j);
% % %     fourier_ratio_acc(j) = acc_dftr_adj(j);
% % % 
%     % instead just assign them as non-adjusted
%     vel_dftr_adj(j) = fourier_ratio_vel(j);
%     acc_dftr_adj(j) = fourier_ratio_acc(j);
%     vel_pct_adj(j) = vel_pct(j);
%     acc_pct_adj(j) = acc_pct(j);
%     
%     %calculate p-values
% %     [dft_p(j) dft_p_vel(j) dft_p_acc(j)] = DFT_perm(gauss_time{j}, gauss_dat{j}, fourier_ratio(j), fourier_ratio_vel(j), fourier_ratio_acc(j), 1000, dftcutoff);
%       dft_p(j) = 999; dft_p_vel(j) = 999; dft_p_acc(j) = 999;
% 
% end
% % end frequency analysis-------------------------------------------------------------------------------------
% % -----------------------------------------------------------------------------------------------------------
% % -----------------------------------------------------------------------------------------------------------
% % -----------------------------------------------------------------------------------------------------------
% % -----------------------------------------------------------------------------------------------------------
% 
% 
% 
% clear data spike_data_withnull temp_spike_data spike_data
% save(['Z:\Users\Yong\FA_for_Chris_' FILE '.mat']);
return;