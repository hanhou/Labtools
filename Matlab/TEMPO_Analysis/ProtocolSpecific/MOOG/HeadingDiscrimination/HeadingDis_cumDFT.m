%-----------------------------------------------------------------------------------------------------------------------
%-- psychometric and neurometric function for heading discrimination task
%--	07/16/04 GY
%-----------------------------------------------------------------------------------------------------------------------

function HeadingDis_cumDFT(data, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);

TEMPO_Defs;
Path_Defs;
ProtocolDefs; 

temp_azimuth = data.moog_params(AZIMUTH,:,MOOG);
temp_elevation = data.moog_params(ELEVATION,:,MOOG);
temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG);
temp_heading   = data.moog_params(HEADING, :, MOOG); 
temp_amplitude = data.moog_params(AMPLITUDE,:,MOOG);
temp_num_sigmas = data.moog_params(NUM_SIGMAS,:,MOOG);
temp_motion_coherence = data.moog_params(COHERENCE,:,MOOG);
temp_spike_rates = data.spike_rates(SpikeChan, :); 
temp_total_trials = data.misc_params(OUTCOME, :);
temp_spike_data = data.spike_data(1,:);   % spike rasters
%now, remove trials from direction and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:length(temp_azimuth);		% a vector of trial indices
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );
stim_type = temp_stim_type( select_trials );
heading = temp_heading( select_trials );
amplitude= temp_amplitude( select_trials );
num_sigmas= temp_num_sigmas( select_trials );
motion_coherence = temp_motion_coherence(select_trials);
spike_rates = temp_spike_rates( select_trials);
total_trials = temp_total_trials( select_trials);
unique_stim_type = munique(stim_type');
unique_heading = munique(heading');
unique_amplitude = munique(amplitude');
unique_num_sigmas = munique(num_sigmas');
unique_motion_coherence = munique(motion_coherence');
disc_heading = unique_heading( floor(length(unique_heading)/2)+1 : end );

%%%%%%%% PSTH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% timebin for plot PSTH
timebin=50;
% sample frequency depends on test duration
frequency=length(temp_spike_data)/length(select_trials);  
% length of x-axis
x_length = frequency/timebin;
% x-axis for plot PSTH
x_time=1:(frequency/timebin);

% remove null trials, bad trials, and trials outside Begtrial~Engtrial
stim_duration = length(temp_spike_data)/length(temp_azimuth);
Discard_trials = find(trials <BegTrial | trials >EndTrial);
for i = 1 : length(Discard_trials)
    temp_spike_data( 1, ((Discard_trials(i)-1)*stim_duration+1) :  Discard_trials(i)*stim_duration ) = 9999;
end
spike_data = temp_spike_data( temp_spike_data~=9999 );
spike_data( find(spike_data>100) ) = 1; % something is absolutely wrong

% count spikes from raster data (spike_data)
max_count = 1;
time_step=1;
for k=1: length(unique_stim_type)
    for i=1:length(unique_heading)
        select = logical( (heading==unique_heading(i)) & (stim_type==unique_stim_type(k)) );            
        act_found = find( select==1 );
        % count spikes per timebin on every same condition trials
        temp_count=[];
        for repeat=1:length(act_found) 
            for n=1:(x_length)
                temp_count(repeat,n)=sum(spike_data(1,(frequency*(act_found(repeat)-1)+time_step):(frequency*(act_found(repeat)-1)+n*timebin)));
                time_step=time_step+timebin;
            end
            time_step=1;                    
        end
        count_y_trial{i,k}(:,:) = temp_count;  % each trial's PSTH 
        dim{i,k} = size(temp_count);
        count_y{i,k} = mean(temp_count);  % take the mean
        max_count_y(i,k) = max(count_y{i,k});
        resp_mean(i,k) = sum(count_y{i,k}(26:58)); % mean firing rate
    end 
    
%     for i=1:length(unique_heading)  
%         % normalize
%   %      count_y{i,k} = count_y{i,k} / max(max_count_y(:,k));
%         
%         % running mean
%         binsize = 3; % take mean of 3 neighbour bins
%         slidesize = 1; % with a slide of 1 bin
%         nn = 1;
%         start_bin = 1;
%         end_bin = start_bin + binsize -1;
%         while end_bin < 100 % the length of count_y{i,k}
%             if sum( count_y{i,k}(start_bin : end_bin) )~=0 % avoid divided by 0 later          
%                 count_y_runmean{i,k}(nn) = mean( count_y{i,k}(start_bin : end_bin) );
%             else 
%                 count_y_runmean{i,k}(nn) = 0;
%             end
%             start_bin = start_bin + slidesize;
%             end_bin = start_bin + binsize -1;
%             nn = nn + 1;
%         end
%     end  
end
  
% now calculate DFT based on responses in each condition (direction and stim_type)
for k = 1 : length(unique_stim_type)
    % find maximum direction
    pref = find( resp_mean(:,k) == max(resp_mean(:,k)) );
    pre(k) = pref(1); % in case two same peaks
%     for i = pre : pre        
%         gdat = count_y_runmean{i,k};
%         dum_dat = gdat;
%         %find hilbert transform to get peak.
%         span = 0;
%         shift_win = 0;
%         peak_time = 40;
%         
%         dummy_dat = gdat(18:63) - mean(gdat(11:20));%only consider 1.5 to 3 sec for finding the peak of the response
%         hil = hilbert(dummy_dat);
%         absabs = abs(hil);
%         HT = abs(hil) - mean(absabs(4:5));% remove the mean to analyze just the peak or else the center of mass will be biased
%         HT(HT<0) = 0;    
%         
%         shift_win_vel = find(hil==max(hil));
%         shift_win_veloc = shift_win_vel(1);
%                 
%         H_time = 1:.05:3;
%         for hh = 1:length(HT)
%             CMsum(hh) = HT(hh)*H_time(hh);
%         end
%         if sum(HT) > 0
%             center_mass = sum(CMsum)/sum(HT);
%             peak_t = find( abs(H_time - center_mass) == min(abs(H_time - center_mass)) );% point which is closest to the center of mass is used.
%         end    
%         
% %         
% %         % stim peak time = 2 sec = 40
% %         shift_win = peak_time - 40; %shift window for fourier transform to get rid of phase shift caused by delayed response.        
%         
% %         %take data from only middle 2 sec (stimulus duration) exlcuding first
% %         %150ms to exclude reation to onset of visual stimulus
% %         for kk=1:40
% %             gaussdata2(kk)=total_data2(kk+20+shift_win); %window is shifted depending on how delayed the peak of the envelope (from the hilbert trans) is
% %         end
% %         gauss_dat2 = gaussdata2;
% %         gaustime=[x_time(21+shift_win):x_time(60+shift_win)];
% %         gausstime=0.05*gaustime; 
% %         %calculate DFT ratio
% %         [f2, amp2, resp_phase2] = FT(gausstime, gauss_dat2, 40, 1, 0);
% %         f1 = mean(amp2(2:4));
% %         f2 = mean(amp2(5:end));
% %         if f2 == 0
% %             f_rat2 = 0;
% %         else
% %             f_rat2 = f1/f2;
% %         end
% %         max_phase2 = mean(resp_phase2(find(amp2 == max(amp2(2:4)))));
% %         mean_phase(i,k) = max_phase2;
% %         fourier_ratio_mean(i,k) = f_rat2;
% %         for j = 1 : dim{i,k}(1)   
% %             total_data = count_y_trial{i,k}(j,1:span);
% %             %if the window is being shifted to extend into the region where the
% %             %recording has stopped, then we will randomly sample from the front of
% %             %the response and fill out the window.
% %             signal = total_data(1:16);
% %             ran_in = randperm(length(signal));
% %             rand_sig = signal(ran_in);
% %             span_new = span;
% % 
% %             if shift_win > 7 & yong_cell % if we shift by more than 7 bins, data cannot be used since there are saccade signials
% %                 if shift_win > 24
% %                     shift_win = 24; %this is as far as we can go for yong's cells
% %                 end
% %                 for ii = 1:(shift_win - 7)
% %                     total_data(span+ii) = rand_sig(ii);
% %                     span_new = span+ii;
% %                 end
% %             end
% % 
% %             %take data from only middle 2 sec (stimulus duration) exlcuding first
% %             %150ms to exclude reation to onset of visual stimulus
% %             for kk=1:40
% %                 gaussdata(kk)=total_data(kk+20+shift_win); %window is shifted depending on how delayed the peak of the envelope (from the hilbert trans) is
% %             end
% %             gauss_dat = gaussdata;
% %             gaustime=[x_time(21+shift_win):x_time(60+shift_win)];
% %             gausstime=0.05*gaustime; 
% %             %calculate DFT ratio
% %             [f, amp, resp_phase] = FT(gausstime, gauss_dat, 40, 1, 0);
% %             f1 = mean(amp(2:4));
% %             f2 = mean(amp(5:end));
% %             if f2 == 0
% %                 f_rat = 0;
% %             else
% %                 f_rat = f1/f2;
% %             end
% %             max_phase = mean(resp_phase(find(amp == max(amp(2:4)))));
% %             fourier_ratio_trial{i,k}(j) = f_rat;  
% %             fourier_ratio_trial_v{i,k}(j) = -f_rat*cos(max_phase2);  % velocity component
% %             fourier_ratio_trial_a{i,k}(j) = f_rat*sin(max_phase2);  % acceleration component
% %             fourier_phase_trial{i,k}(j)=max_phase*180/pi;
% % %             if max_phase2*180/pi>-90 & max_phase2*180/pi<90
% % %                fourier_ratio_trial_sign{i,k}(j) = -f_rat; 
% % %             else
% % %                fourier_ratio_trial_sign{i,k}(j) = f_rat;
% % %             end            
% %          end
%     end
 end

%     %%%%%%% replace spike_rates with fourier_ratio_trial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % save a copy for spike_rates first
% for nn = 1:2    
%     spike_rates_copy = spike_rates;
%     for k = 1 : length(unique_stim_type)
%         for i = 1 : length(unique_heading)
%             select = logical(heading==unique_heading(i) & stim_type==unique_stim_type(k));
%             act_find = find(select==1);
%             for j = 1 : length(act_find)
%                 if nn==1
%                    spike_rates(act_find(j)) = fourier_ratio_trial_v{i,k}(j); 
%                 else 
%                    spike_rates(act_find(j)) = fourier_ratio_trial_a{i,k}(j);
%                 end
%             end
%         end
%     end
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %neurometric dataset and calculate ROC, Choice Probability(CP)
%     %determine for each trial whether monkey chooses leftward(target1) or rightward(tarket2)    
%     LEFT = 1;
%     RIGHT = 2;
%     for i= 1 : length(spike_rates) 
%         temp = data.event_data(1,:,i + BegTrial-1);
%         events = temp(temp>0);  % all non-zero entries
%         if (sum(events == IN_T1_WIN_CD) > 0)
%             choice(i) = RIGHT;
%         elseif (sum(events == IN_T2_WIN_CD) > 0)
%             choice(i) = LEFT;
%         else
%          %   choice(i) = RIGHT;
%             disp('Neither T1 or T2 chosen.  This should not happen!.  File must be bogus.');
%         end
%     end
% 
%     one_repetition = length(unique_heading)*length(unique_stim_type);
%     repetition = floor( length(spike_rates)/one_repetition ); % take minimum repetition
%     resp_heading = [];
%     Z_Spikes = spike_rates;
%     % z-score data for later cp analysis across headings
%     for k = 1:length(unique_stim_type)
%         for i = 1:length(unique_heading)
%             select =logical( (heading == unique_heading(i)) & (stim_type == unique_stim_type(k)) ) ;  
%             z_dist = spike_rates(select);
%             z_dist = (z_dist - mean(z_dist))/std(z_dist);
%             Z_Spikes(select) = z_dist;
%         end
%     end
%     Z_Spikes_Ori = Z_Spikes; % keep a Z_Spikes unchanged for later use
%     % now group neuronal data into two groups according to monkey's choice
%     for k = 1:length(unique_stim_type)    % notice that the condition is double than disc_heading    
%         for i = 1:length(unique_heading)
%             select =logical( (heading == unique_heading(i)) & (stim_type == unique_stim_type(k)) ) ;  
%             resp{k,i} = spike_rates(select);
% 
%             resp_mat_copy{k}(i) = mean(spike_rates_copy(select)); % metric based on mean firing rate
%             resp_mat{k}(i) = mean(resp{k,i});     % use the DFT for each trial     
%             resp_mat_err{k}(i) = std(resp{k,i}) / sqrt(repetition);
% 
%             % calculate CP, group data based on monkey's choice 
%             resp_left_choose{k,i} = spike_rates(select & (choice == LEFT) );
%             resp_right_choose{k,i} = spike_rates(select & (choice == RIGHT) );
%             if (length(resp_left_choose{k,i}) <= 3) | (length(resp_right_choose{k,i}) <= 3)   % make sure each condition has at least 3 data values
%                 Z_Spikes(select) = 9999;   % similar to NaN, just make a mark            
%                 CP{k}(i) = 9999;
%             else
%                 CP{k}(i) = 0;
%             end 
%         end  
%         % now cross all data 
%         resp_left_all{k} = Z_Spikes(  (stim_type == unique_stim_type(k)) & (choice == LEFT) & (Z_Spikes~=9999) ); 
%         resp_right_all{k} = Z_Spikes(  (stim_type == unique_stim_type(k)) & (choice == RIGHT) & (Z_Spikes~=9999) ); 
%         resp_all{k} = Z_Spikes(  (stim_type == unique_stim_type(k)) & (Z_Spikes~=9999) ); 
%     end
% 
%     % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %    % decide whether ves and vis is congruent tuning. Fit line by linear
%     % 	% regression first and compare the sign of each condition to decide whether
%     % 	% congruent or opposite, this is used to check whether congruent cells lead
%     % 	% to better neuronal performance in combined condition, and vice versa
%     for k = 1 : length(unique_stim_type)
%         [rr,pp] = corrcoef(unique_heading, resp_mean(:,k) );  % slope based on firing rate
%         line_re{k} = rr(1,2);
%         line_p{k} = pp(1,2);
%         p{k} = polyfit(unique_heading, resp_mean(:,k),1 );    
% 
%       %  [rr_DFT_trial,pp_DFT_trial] = corrcoef(unique_heading, resp_mat{k}(:));  % slope based on DFT trials
%         [rr_DFT_trial,pp_DFT_trial] = corrcoef(unique_heading, fourier_ratio_mean(:,k) );
%         line_re_DFT_trial{k} = rr_DFT_trial(1,2);
%         line_p_DFT_trial{k} = pp_DFT_trial(1,2);
%         p_DFT_trial{k} = polyfit(unique_heading, resp_mat{k}(:),1 );    
%     end
% 
%     % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % now calculate propotion correct from area under ROC curves, each heading is compared to 0 heading
%     % neurothreshold 
%     fit_data_neuro = [];
%     fit_data_neuro_cut = [];
%     for k = 1 : length(unique_stim_type)
%         for i = 1 : length(unique_heading)-1   % subtract the 0 heading
%             trials_n =logical( (heading == unique_heading(i)) & (stim_type == unique_stim_type(k)) ) ;
%             fit_data_neuro_cum{k}(i,3) = sum(trials_n);  % for later function fit use
%             if i < (1+length(unique_heading))/2
%                  % in rocN, prefer first, but since plot rightward choice, thus reverse, here we asume leftward is the preferred direction            
%                  Neuro_correct{k}(i) =  rocN( resp{k,(1+length(unique_heading))/2},resp{k,i},100 ); % compare to the 0 heading condition, which is straght ahead
%              else
%                  Neuro_correct{k}(i) =  rocN( resp{k,(1+length(unique_heading))/2}, resp{k,(i+1)},100 ); % compare to the 0 heading condition, which is straght ahead
%              end
%              if line_re_DFT_trial{k} > 0 % if left heading is not the larger responses, then linear regression will be positive 
%                  % we don't know whether left heading is for sure smaller than right heading,thus ROC curve might be flipped above or below the unity line
%                  % here we asume if left response larger than right response then asign the left to be preferred direction
%                  Neuro_correct{k}(i) = 1 - Neuro_correct{k}(i);            
%              end  
%          end
%     end
%     % choice probability
%     for k = 1 : length(unique_stim_type)
%         for i = 1 : length(unique_heading)  
%             if CP{k}(i)~=9999
%                CP{k}(i) = rocN( resp_left_choose{k,i},resp_right_choose{k,i},100 );
%             else
%                CP{k}(i) = NaN;
%             end
%             if (length(resp_left_all{k}) > 3) | (length(resp_right_all{k}) > 3)
%             %if  (length(resp_left_all{k}) / length(resp_all{k}) >0.25) |  (length(resp_left_all{k}) / length(resp_all{k}) <0.75)
%                 CP_all{k} = rocN( resp_left_all{k},resp_right_all{k},100 );
%             else
%                 CP_all{k} = NaN; 
%             end
%             if  line_re_DFT_trial{k} > 0
%                 CP{k}(i) = 1 - CP{k}(i);
%                 CP_all{k} = 1 - CP_all{k};
%             end
%         end
%     end
%   
%     %%%%%% use Wichman's MLE method to estimate threshold and bias
%     % for k = 1:length(unique_stim_type)
%     for k = 1:1    
%         fit_data_neuro_cum{k}(:,1) = unique_heading(unique_heading~=0);
%         fit_data_neuro_cum{k}(:,2) = Neuro_correct{k}(:);
% %        wichman_neu = pfit(fit_data_neuro_cum{k},'plot_opt','no plot','shape','cumulative gaussian','n_intervals',1,'FIX_LAMBDA',0.001,'sens',0,'compute_stats','false','verbose','false');
% %        Thresh_neu{k} = wichman_neu.params.est(2);
%         if sum(Neuro_correct{k}(:)) < 0.1 | sum(Neuro_correct{k}(:)) >7.9
%             Thresh_neu{k} = 888;
%         else
%             [bb,tt] = cum_gaussfit_max1(fit_data_neuro_cum{k});
%             Thresh_neu{k} = tt;
%             if Thresh_neu{k}<0 | Thresh_neu{k}> 300
%                 Thresh_neu{k} = 300;           
%             end
%         end
%     end
%     Neu(nn) = Thresh_neu{1};
%     CPP(nn) = CP_all{1};
%     Rocvalue{nn} = Neuro_correct{1};
% end

%% plot PSTH now %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get the largest count_y so that make the scale in each figures equal    
% plot two lines as stimulus start and stop marker
x_start = [StartEventBin(1,1)/timebin, StartEventBin(1,1)/timebin];
x_stop =  [StopEventBin(1,1)/timebin,  StopEventBin(1,1)/timebin];
y_marker=[0,max(max(max_count_y))];
% define figure
figure(2);
set(2,'Position', [5,5 1000,680], 'Name', 'Tuning');
orient landscape;
axis off;

xoffset=0;
yoffset=0;

for k=1: length(unique_stim_type)     
    for i=1:length(unique_heading)        
        axes('position',[0.31*(k-1)+0.1 (0.92-0.09*i) 0.25 0.08]);          
        bar( x_time,count_y{i,k}(1,:) );  
        hold on;
        plot( x_start, y_marker, 'r-');
        plot( x_stop,  y_marker, 'r-');
        set( gca, 'xticklabel', ' ' );
        % set the same scale for all plot
        xlim([0,x_length]);
        ylim([0,max(max(max_count_y))]);
    end 
    axes('position',[0 0 1 1]); 
    xlim([-50,50]);
    ylim([-50,50]);
    text(-25,45, 'vestibular'); 
    text(0,45, 'visual'); 
    text(25,45, 'combined'); 
    text(-45, 45, FILE);
%     for j = 1:length(unique_heading)
%         text(-50+k*30, 48-j*9, num2str(resp_mean(j,k)) );
%         text(-50+k*30, 46-j*9, num2str(fourier_ratio(j,k)) );
%         text(-50+k*30, 44-j*9, num2str(fourier_phase(j,k)) );
%         if k == 1 
%             text(-48, 45-j*9, num2str(unique_heading(j)) );
%         end
%     end 
    axis off;
    hold on;
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Also, write out some summary data to a cumulative summary file
sprint_txt = ['%s']; 
for i = 1 : 300 % this should be large enough to cover all the data that need to be exported
     sprint_txt = [sprint_txt, ' %4.3f'];    
end
%buff = sprintf(sprint_txt, FILE, shift_win_veloc, count_y_runmean{pre,1}, abs(hil)  );
buff = sprintf(sprint_txt, FILE, count_y{pre(1),1},count_y{pre(2),2},count_y{pre(3),3}  );
%buff = sprintf(sprint_txt, FILE, p{1}, p{2}, p{3},line_re{:},line_p{:}, p_DFT_MEAN{1}, line_p_DFT_MEAN{1},p_DFT_trial{1},line_p_DFT_trial{1} );
outfile = [BASE_PATH 'ProtocolSpecific\MOOG\HeadingDiscrimination\HeadingDiscri_cumPSTH_CueCombine.dat'];
printflag = 0;
if (exist(outfile, 'file') == 0)   % file does not yet exist
    printflag = 1;
end
fid = fopen(outfile, 'a');
if (printflag)
    fprintf(fid, 'FILE\t         theshold\t CP\t Rfiring\t Rdfttrial\t Rdftmean\t phase\t ratio\t');
    fprintf(fid, '\r\n');
end
fprintf(fid, '%s', buff);
fprintf(fid, '\r\n');
fclose(fid);
%---------------------------------------------------------------------------------------
return;