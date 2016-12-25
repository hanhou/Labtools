% analysis of eyetrace for heading discrimination
% including eye position, eye velocity
% -GY
function HeadingDis_cum_eyetrace(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

TEMPO_Defs;
Path_Defs;
ProtocolDefs; %contains protocol specific keywords - 1/4/01 YG

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
trials = 1:length(temp_azimuth);
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );
azimuth = temp_azimuth(select_trials);
elevation = temp_elevation(select_trials);
stim_type = temp_stim_type(select_trials);
amplitude = temp_amplitude(select_trials);
heading = temp_heading( select_trials );
spike_rates = temp_spike_rates(select_trials);
eye_data = data.eye_data(:,:,select_trials);

%spike_rates(1:length(spike_rates)) = rand(1, length(spike_rates)); % for staircase and psychophysical

unique_azimuth = munique(azimuth');
unique_elevation = munique(elevation');
unique_stim_type = munique(stim_type');
unique_amplitude = munique(amplitude');
unique_heading = munique(heading');

one_repetition = length(unique_heading)*length(unique_stim_type);
repetition = floor( length(spike_rates)/one_repetition ); % take minimum repetition


%% Plot eyetrace. HH20150521

eyeX = squeeze(data.eye_data(1,:,:));
eyeY = squeeze(data.eye_data(2,:,:));

trial_begin = 1000/5;
trial_end = 3000/5;

eyeX = eyeX(trial_begin:trial_end,:);
eyeY = eyeY(trial_begin:trial_end,:);

figure();
plot(eyeX(:),eyeY(:));





%%

% for i = 1:length(spike_rates)
%     eyeplot(i,:) = eye_data(1,:,i);
% end
% figure;
% plot(eyeplot','-');
% xlim([201,700]);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % calculate eye velocity based on trials to replace spike_rate and then
% % calculate neuronal threshold and CP, similar to the analysis of accelerometer
% if sum( eye_data(1,201:300,1) ) ~= 0    % take the left coil if it exist
%     eye_chan = 1;
% else
%     eye_chan = 3;
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %neurometric dataset and calculate ROC, Choice Probability(CP)
% %determine for each trial whether monkey chooses leftward(target1) or rightward(tarket2)    
% LEFT = 1;
% RIGHT = 2;
% for i= 1 : length(spike_rates) 
%     temp = data.event_data(1,:,i + BegTrial-1);
%     events = temp(temp>0);  % all non-zero entries
%     if (sum(events == IN_T1_WIN_CD) > 0)
%         choice(i) = RIGHT;
%     elseif (sum(events == IN_T2_WIN_CD) > 0)
%         choice(i) = LEFT;
%     else
%         disp('Neither T1 or T2 chosen.  This should not happen!.  File must be bogus.');
%     end
% end
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% raw eye traces
% % % stimulus
% % for i = 1 : length(unique_heading )
% %     select = find( heading==unique_heading(i) & stim_type==1 ); % only vestibular condition
% %     for n = 1 : length(select)  
% % %             % for vergence
% % %             DC_begin = mean( (eye_data(1,201:300,select(n))-eye_data(3,201:300,select(n))) );
% % %             DC_end = mean( (eye_data(1,501:600,select(n))-eye_data(3,501:600,select(n))) );
% % %             DC = mean([DC_begin, DC_end]);            
% % %             verg_azi(n) = mean( (eye_data(1,301:500,select(n))-eye_data(3,301:500,select(n))) )- DC;
% % %             verg_trace_trial(n,:) = eye_data(1,201:600,select(n))-eye_data(3,201:600,select(n)) - DC;
% % %             verg_trial{i,k}(n) = mean(verg_trace_trial(n,51:250)); % the mean middle 1 sec
% %         DC_hor = mean( eye_data(eye_chan,201:300,select(n)) );
% %         eye_pos_hor_trial(n,:) = eye_data(eye_chan,201:700,select(n)) - DC_hor;     
% %         DC_ver = mean( eye_data(eye_chan+1,201:300,select(n)) );
% %         eye_pos_ver_trial(n,:) = eye_data(eye_chan+1,201:700,select(n)) - DC_ver; 
% %     end
% %     eye_pos_hor_head(i,:) = median( eye_pos_hor_trial(:,:) );
% %     eye_pos_ver_head(i,:) = median( eye_pos_ver_trial(:,:) );        
% % end
% % % now grouped by choice 
% % for i = 1: round(length(unique_heading)/2) % only vestibular condition
% %     select_left = find( stim_type==1 & choice==LEFT & (heading==unique_heading(i) | heading==unique_heading(length(unique_heading)+1-i)) );
% %     select_right = find( stim_type==1 & choice==RIGHT & (heading==unique_heading(i) | heading==unique_heading(length(unique_heading)+1-i)) );
% %     if length(select_left) >= 1;
% %         for n = 1: length(select_left)
% %             DC_hor_left = mean( eye_data(eye_chan,201:300,select_left(n)) );
% %             eye_hor_left_trial(n,:) = eye_data(eye_chan,201:700,select_left(n)) - DC_hor_left;            
% %         end
% %         eye_pos_hor_left(i,:) = median( eye_hor_left_trial(:,:) );        
% %     else
% %         eye_pos_hor_left(i,1:500) = 99; % mark         
% %     end
% %     
% %     if length(select_right) >= 1;
% %         for n = 1: length(select_right)
% %             DC_hor_right = mean( eye_data(eye_chan,201:300,select_right(n)) );
% %             eye_hor_right_trial(n,:) = eye_data(eye_chan,201:700,select_right(n)) - DC_hor_right;            
% %         end
% %         eye_pos_hor_right(i,:) = median( eye_hor_right_trial(:,:) );        
% %     else
% %         eye_pos_hor_right(i,1:500) = 99; % mark         
% %     end
% % end
% % 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for win_shift = 6:6 % shift window by every 20 bins, which is 100 ms, 6th should be the middle one sec
%   
%     for i = 1 : length(spike_rates)
%         DC = mean( eye_data(1,201:300,i)-eye_data(3,201:300,i) );
%         trace_trial = eye_data(eye_chan,201:700,i) - mean(eye_data(eye_chan,201:300,i)); % normalize to start about 0
%         trace_vergence = eye_data(1,201:600,i) - eye_data(3,201:600,i) -DC;
%         vel_trial = diff(trace_trial)/0.005; % velocity trace
%         vel_trial(500) = vel_trial(499);
%         resp_eyevelocity(i) = sum(vel_trial(1+20*(win_shift-1):200+20*(win_shift-1)));
%    %     resp_eyeposition(i) = sum( trace_trial(1+20*(win_shift-1):200+20*(win_shift-1)) );  
%         resp_eyeposition(i) = sum( trace_vergence(1+20*(win_shift-1):200+20*(win_shift-1)) );  
%     end
%     % z-score data for later cp analysis across headings
%     for k = 1:length(unique_stim_type)
%         for i = 1:length(unique_heading)
%             select =logical( (heading == unique_heading(i)) & (stim_type == unique_stim_type(k)) ) ;  
%             % for firing rate
%             z_dist = spike_rates(select);
%             z_dist = (z_dist - mean(z_dist))/std(z_dist);
%             Z_Spikes(select) = z_dist;
% 
%             % for eye velocity
%             z_dist_vel = resp_eyevelocity(select);
%             z_dist_vel = (z_dist_vel - mean(z_dist_vel))/std(z_dist_vel);
%             Z_Spikes_vel(select) = z_dist_vel;    
%             
%             % for eye position
%             z_dist_pos = resp_eyeposition(select);
%             z_dist_pos = (z_dist_pos - mean(z_dist_pos))/std(z_dist_pos);
%             Z_Spikes_pos(select) = z_dist_pos;    
%         end
%     end
% 
%     % now group neuronal data into two groups according to monkey's choice
%     for k = 1:1    % notice that the condition is double than disc_heading    
%         resp_co = [];
%         resp_co_vel = [];
%         resp_co_pos = [];
%         for i = 1:length(unique_heading)
%             select =logical( (heading == unique_heading(i)) & (stim_type == unique_stim_type(k)) ) ;  
%             resp{k,i} = spike_rates(select); 
%             resp_co = [resp_co, resp{1,i}];
%             resp_mat{k}(i) = mean(resp{k,i});          
%             resp_mat_err{k}(i) = std(resp{k,i}) / sqrt(repetition);
% 
%             resp_eyevel{k,i} = resp_eyevelocity(select);    
%             resp_co_vel = [resp_co_vel, resp_eyevel{1,i}];
%             resp_mat_vel{k}(i) = mean(resp_eyevel{k,i});          
%             resp_mat_vel_err{k}(i) = std(resp_eyevel{k,i}) / sqrt(repetition);
%             
%             resp_eyepos{k,i} = resp_eyeposition(select);    
%             resp_co_pos = [resp_co_pos, resp_eyepos{1,i}];
%             resp_mat_pos{k}(i) = mean(resp_eyepos{k,i});          
%             resp_mat_pos_err{k}(i) = std(resp_eyepos{k,i}) / sqrt(repetition);
% 
%             % calculate CP, group data based on monkey's choice 
% %             resp_left_choose_vel{k,i} = resp_eyevelocity(select & (choice == LEFT) );  % eye vel
% %             resp_right_choose_vel{k,i} = resp_eyevelocity(select & (choice == RIGHT) );  % eye vel
% 
%             resp_left_choose{k,i} = spike_rates(select & (choice == LEFT) );    % spike_rates
%             resp_right_choose{k,i} = spike_rates(select & (choice == RIGHT) );    %spike_rates
% 
%             if (length(resp_left_choose{k,i}) <= 3) | (length(resp_right_choose{k,i}) <= 3)   % make sure each condition has at least 3 data values
%                 Z_Spikes_vel(select) = 9999;   % similar to NaN, just make a mark 
%                 Z_Spikes_pos(select) = 9999; 
%                 Z_Spikes(select) = 9999;   
%             end 
%         end  
%         % now cross all data 
%         select_left = (stim_type == unique_stim_type(k)) & (choice == LEFT) & (Z_Spikes~=9999) & (Z_Spikes_vel~=9999) & (Z_Spikes_pos~=9999);
%         select_right = (stim_type == unique_stim_type(k)) & (choice == RIGHT) & (Z_Spikes~=9999)& (Z_Spikes_vel~=9999) & (Z_Spikes_pos~=9999);
%         resp_left_all_vel{k} = Z_Spikes_vel(  select_left );   % eye vel
%         resp_right_all_vel{k} = Z_Spikes_vel(  select_right );  % eye vel 
%         
%         resp_left_all_pos{k} = Z_Spikes_pos(  select_left );   % eye pos
%         resp_right_all_pos{k} = Z_Spikes_pos(  select_right );  % eye pos
%         
%         resp_left_all{k} = Z_Spikes(  select_left );     % spike rate
%         resp_right_all{k} = Z_Spikes(  select_right );    % spike rate
% 
%         % calculate slope
%         [rr,pp] = corrcoef(unique_heading, resp_mat_vel{k}(:));
%         line_re_vel{k} = rr(1,2);
%         line_p_vel{k} = pp(1,2);    
%         
%         [rr,pp] = corrcoef(unique_heading, resp_mat_pos{k}(:));
%         line_re_pos{k} = rr(1,2);
%         line_p_pos{k} = pp(1,2);  
% 
%         [rr,pp] = corrcoef(unique_heading, resp_mat{k}(:));
%         line_re{k} = rr(1,2);
%         line_p{k} = pp(1,2);     
%         
%         % calculate and remove correlation between neural activity and eye velocity        
%         s = polyfit(resp_co_vel, resp_co,1);      
%         for i = 1:length(unique_heading)
%             resp_vel_correct{k,i} = resp{k,i}- resp_eyevel{k,i}*s(1)+s(2);            
%         end
%         % calculate and remove correlation between neural activity and eye
%         % position
%         s = polyfit(resp_co_pos, resp_co,1);      
%         for i = 1:length(unique_heading)
%             resp_pos_correct{k,i} = resp{k,i}- resp_eyepos{k,i}*s(1)+s(2);            
%         end
%     end
% 
%     %-------------------------------------------------------------------------
%     % calculate threshold
%     for k = 1 : 1
%         for i = 1 : length(unique_heading)-1   % subtract the 0 heading
%             trials_n =logical( (heading == unique_heading(i)) & (stim_type == unique_stim_type(k)) ) ;
%             fit_data_neuro_cum{k}(i,3) = sum(trials_n);  % for later function fit use
%             if i < round(length(unique_heading)/2)
%                  Neuro_correct{k}(i) =  rocN( resp_pos_correct{k,length(unique_heading)-i+1},resp_pos_correct{k,i},100 );    % antineuron model
%              else
%                  Neuro_correct{k}(i) =  rocN( resp_pos_correct{k,length(unique_heading)-i}, resp_pos_correct{k,(i+1)},100 );   % antineuron model  
%              end
%              if line_re{k} > 0
%                  Neuro_correct{k}(i) = 1 - Neuro_correct{k}(i);            
%              end  
%          end
%     end
%     % choice probability
%     for k = 1 : 1  
%         if (length(resp_left_all{k}) > 3) & (length(resp_right_all{k}) > 3)
%             CP_vel{k} = rocN( resp_right_all_vel{k},resp_left_all_vel{k},100 ); % right is positve
%             CP_pos{k} = rocN( resp_right_all_pos{k},resp_left_all_pos{k},100 ); % right is positve            
%             CP_all{k} = rocN( resp_left_all{k},resp_right_all{k},100 );
%         else
%             CP_vel{k} = NaN; 
%             CP_pos{k} = NaN;
%             CP_all{k} = NaN;
%         end
% 
%         if line_re{k} > 0
%             CP_all{k} = 1 - CP_all{k};               
%         end 
%     end    
%     %%%%%% use Wichman's MLE method to estimate threshold and bias
%     for k = 1:1   % calculate vestibular temporarily
%         fit_data_neuro_cum{k}(:,1) = unique_heading(unique_heading~=0);
%         fit_data_neuro_cum{k}(:,2) = Neuro_correct{k}(:);
%         wichman_neu = pfit(fit_data_neuro_cum{k},'plot_opt','no plot','shape','cumulative gaussian','n_intervals',1,'FIX_LAMBDA',0.001,'sens',0,'compute_stats','false','verbose','false');
%         Thresh_neu{k} = wichman_neu.params.est(2);
%         Bias_neu{k} = wichman_neu.params.est(1);
%         neu_perf{k} = [wichman_neu.params.est(1),wichman_neu.params.est(2)];
%         if Thresh_neu{k}>300 | Thresh_neu{k}<0
%             Thresh_neu{k}=300;
%         end
%     end
%     CP_vel_shift(win_shift) = CP_vel{1};  
%     CP_pos_shift(win_shift) = CP_pos{1};
%     Thresh_neu_shift(win_shift) = Thresh_neu{1};
% end
% % % %--------------------------------------------------------------------------
% % %do permutation to test the significance of CP_all{k}, re-calculate CP 1000 times
% % perm_num = 1000;
% % Z_Spikes_perm = Z_Spikes_vel;
% % bin = 0.005;
% % x_bin = 0 : bin : 1;
% % for k = 1:1 % vestibular temporarily    
% %     for n = 1 : perm_num
% %         % temperarilly only use near-threshold heading angles where monkey make a guess mainly
% %         select = logical( (stim_type == unique_stim_type(k)) & (Z_Spikes~=9999) );
% %         Z_Spikes_con{k} = Z_Spikes_perm( select );
% %         Z_Spikes_con{k} = Z_Spikes_con{k}(randperm(length(Z_Spikes_con{k})));   % permute spike_rates
% %         Z_Spikes_perm(select) = Z_Spikes_con{k};    % now in spike_rates, the corresponding data were permuted already
% %         
% %         resp_left_all_perm{k} = Z_Spikes_perm(  (stim_type == unique_stim_type(k)) & (choice == LEFT) & (Z_Spikes~=9999) ); 
% %         resp_right_all_perm{k} = Z_Spikes_perm(  (stim_type == unique_stim_type(k)) & (choice == RIGHT) & (Z_Spikes~=9999) ); 
% %     
% %         if  (length(resp_left_all{k}) > 3) & (length(resp_right_all{k}) > 3) 
% %             CP_all_perm{k}(n) = rocN( resp_left_all_perm{k}, resp_right_all_perm{k},100 );
% %         else
% %             CP_all_perm{k}(n) = NaN; 
% %         end
% %         if  line_re_vel{k} > 0       
% %             CP_all_perm{k}(n) = 1 - CP_all_perm{k}(n);             
% %         end  
% %     end
% %     % now calculate p value or significant test
% %     if (length(resp_left_all{k}) > 3) & (length(resp_right_all{k}) > 3) 
% %         hist_perm(k,:) = hist( CP_all_perm{k}(:), x_bin );  % for permutation
% %         bin_sum = 0;
% %         n = 0;
% %         while ( n < (CP_vel{k}/bin) )
% %              n = n+1;
% %              bin_sum = bin_sum + hist_perm(k, n);
% %              if CP_vel{k} > 0.5                  % note it's two tail test
% %                 p{k} = 2*(perm_num - bin_sum)/ perm_num;    % calculate p value for CP_all
% %              else
% %                 p{k} = 2* bin_sum / perm_num;
% %              end
% %         end
% %     else
% %         p{k} = NaN;
% %     end 
% % end
% % 
%  
% %%% Ancova analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % vestibular temporally
% xx(1:length(resp_left_all_pos{1}),1) = resp_left_all_pos{1};
% xx(length(resp_left_all_pos{1})+1:length(resp_left_all_pos{1})+length(resp_right_all_pos{1}),1) = resp_right_all_pos{1};
% yy(1:length(resp_left_all{1}),1) = resp_left_all{1}';
% yy(length(resp_left_all{1})+1:length(resp_left_all{1})+length(resp_right_all{1}),1) = resp_right_all{1}';
% cc(1:length(resp_left_all{1}),1) = 1; % mark left or right
% cc(length(resp_left_all{1})+1:length(resp_left_all{1})+length(resp_right_all{1}),1) = 2; % mark left or right
% [h1,atab1,ctab1,stats1]=aoctool(xx,yy,cc,'','','','','off'); % fit seperate lines
% [h2,atab2,ctab2,stats2]=aoctool(xx,yy,cc,'','','','','off','same line'); % fit same line to calculate residures
% % remove correlation between eye signal and neural activity, calculate resudules to re-calculate CP
% resp_yy = yy - (ctab2{3,2}*xx + ctab2{2,2});
% resp_left_all_yy = resp_yy(1:length(resp_left_all{1}));
% resp_right_all_yy = resp_yy(length(resp_left_all{1})+1:length(resp_left_all{1})+length(resp_right_all{1}));
% % recalculate CP
% if (length(resp_left_all{1}) > 3) & (length(resp_right_all{1}) > 3)
%     CP_residule = rocN( resp_left_all_yy,resp_right_all_yy, 100 );
% else
%     CP_residule = NaN;
% end
% if line_re{1} > 0
%     CP_residule = 1 - CP_residule;
% end
% % 
% % figure(3);
% % plot(xx(1:length(resp_left_all{1}),1),yy(1:length(resp_left_all{1}),1),'bo');
% % hold on;
% % plot(xx(length(resp_left_all{1})+1:length(resp_left_all{1})+length(resp_right_all{1}),1),yy(length(resp_left_all{1})+1:length(resp_left_all{1})+length(resp_right_all{1}),1),'ro');
% % plot([0,0],[min(yy), max(yy)],'k--');
% % text(min(xx),max(yy),num2str(atab1{2,6}) );
% % text(min(xx),max(yy)-0.5,num2str(atab1{3,6}) );
% % text(min(xx),max(yy)-1,num2str(atab1{4,6}) );
% % text(min(xx),max(yy)-2,num2str(CP_all{1}) );
% % text(min(xx),max(yy)-2.5,num2str(CP_residule) );
% % text(min(xx),max(yy)-3,num2str(CP_vel{1}) );
% % xlim([min(xx)-2,max(xx)]);
% % ylim([min(yy), max(yy)]);
% % xlabel('eye velocity (deg/s)');
% % ylabel('Z-scored neural activity');
% % title(FILE);
% 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %output to text file
% sprint_txt = ['%s'];
% for i = 1 : 2000
%      sprint_txt = [sprint_txt, ' %1.3f'];    
% end
% % 
% % buff = sprintf(sprint_txt, FILE, eye_pos_hor_head(1,:),eye_pos_hor_head(2,:),eye_pos_hor_head(3,:),...
% %                                            eye_pos_hor_head(4,:),eye_pos_hor_head(5,:),eye_pos_hor_head(6,:),...
% %                                            eye_pos_hor_head(7,:),eye_pos_hor_head(8,:),eye_pos_hor_head(9,:),...
% %                                            eye_pos_ver_head(1,:),eye_pos_ver_head(2,:),eye_pos_ver_head(3,:),...
% %                                            eye_pos_ver_head(4,:),eye_pos_ver_head(5,:),eye_pos_ver_head(6,:),...
% %                                            eye_pos_ver_head(7,:),eye_pos_ver_head(8,:),eye_pos_ver_head(9,:),...
% %                                            eye_pos_hor_left(1,:),eye_pos_hor_left(2,:),eye_pos_hor_left(3,:),eye_pos_hor_left(4,:),eye_pos_hor_left(5,:),...
% %                                            eye_pos_hor_right(1,:),eye_pos_hor_right(2,:),eye_pos_hor_right(3,:),eye_pos_hor_right(4,:),eye_pos_hor_right(5,:) );
% % outfile = [BASE_PATH 'ProtocolSpecific\MOOG\HeadingDiscrimination\Eye_velocity_raw.dat'];
% % printflag = 0;
% % if (exist(outfile, 'file') == 0)    %file does not yet exist
% %     printflag = 1;
% % end
% % fid = fopen(outfile, 'a');
% % if (printflag)
% %     fprintf(fid, 'FILE\t');
% %     fprintf(fid, '\r\n');
% % end
% % fprintf(fid, '%s', buff);
% % fprintf(fid, '\r\n');
% % fclose(fid);
% 
% %buff2 = sprintf(sprint_txt, FILE, Thresh_neu{1}, CP_vel{1},p{1},CP_vel_neu{1}, atab1{2,6},atab1{3,6},atab1{4,6}, CP_all{1},CP_residule );
% %buff2 = sprintf(sprint_txt, FILE, Thresh_neu_shift(win_shift), CP_residule, atab1{2,6},atab1{3,6},atab1{4,6});
% buff2 = sprintf(sprint_txt,FILE, mean_accy, mean_accx);
% outfile2 = [BASE_PATH 'ProtocolSpecific\MOOG\HeadingDiscrimination\accelerometerraw.dat'];
% printflag = 0;
% if (exist(outfile2, 'file') == 0)    %file does not yet exist
%     printflag = 1;
% end
% fid = fopen(outfile2, 'a');
% if (printflag)
%     fprintf(fid, 'FILE\t');
%     fprintf(fid, '\r\n');
% end
% fprintf(fid, '%s', buff2);
% fprintf(fid, '\r\n');
% fclose(fid);

return;

