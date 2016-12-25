%-----------------------------------------------------------------------------------------------------------------------
%-- accelerometer measurement for psychometric and neurometric function for heading discrimination task
%--	07/23/04
%-----------------------------------------------------------------------------------------------------------------------

function Accelerometer_cum_1D_Jing(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

TEMPO_Defs;
Path_Defs;
ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP

temp_azimuth = data.moog_params(AZIMUTH,:,MOOG);
temp_elevation = data.moog_params(ELEVATION,:,MOOG);
temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG);
temp_heading   = data.moog_params(HEADING, :, MOOG); 
temp_amplitude = data.moog_params(AMPLITUDE,:,MOOG);
temp_num_sigmas = data.moog_params(NUM_SIGMAS,:,MOOG);
temp_motion_coherence = data.moog_params(COHERENCE,:,MOOG);
temp_spike_rates = data.spike_rates(SpikeChan, :); 
temp_total_trials = data.misc_params(OUTCOME, :);
temp_synpluse = data.spike_data(2,:,:);
temp_accel = data.eye_data(5,:,:);
temp_event = data.event_data(1,:,:);
trials = 1:length(temp_heading);		% a vector of trial indices
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );

% stim_type = temp_stim_type( select_trials );
% heading = temp_heading( select_trials );
% amplitude= temp_amplitude( select_trials );
% num_sigmas= temp_num_sigmas( select_trials );
% motion_coherence = temp_motion_coherence(select_trials);
% spike_rates = temp_spike_rates( select_trials);
total_trials = temp_total_trials( select_trials);
% 
% unique_stim_type = munique(stim_type');
% unique_heading = munique(heading');
% unique_amplitude = munique(amplitude');
% unique_num_sigmas = munique(num_sigmas');
% unique_motion_coherence = munique(motion_coherence');

% one_repetition = length(unique_heading)*length(unique_stim_type);
% repetition = floor( length(spike_rates)/one_repetition ); % take minimum repetition
% replace spike_rate with accelerometer data
%calculate mean offset first, take the first second as the control
offset_y = mean(mean(data.eye_data(5,1:200,:)));
offset_y_foraft = mean(mean(data.eye_data(6,1:200,:)));

% j=1;
% for  i=1:2000
%     if temp_azimuth(i) == 0
%         synpluse(:,j)=temp_synpluse(1,:,i);
%         accel(:,j)=temp_accel(1,:,i)-offset_y;
%         event(:,j)=temp_event(1,:,i);
%         j=j+1;
%     end;    
% end

for i=1:length(total_trials)
    synpluse(:,i)=temp_synpluse(1,:,i);
    accel(:,i)=temp_accel(1,:,i)-offset_y;
    event(:,i)=temp_event(1,:,i);
end

for i=1:length(total_trials)
    for j=1:5000
        if synpluse(j,i)==1
            firstSyn(i)=j;
            break;
        end;
    end;
       
    visualBegintime(i) = find(event(:,i) == 4);
    startInd=round(firstSyn(i)/5);
    acc(:,i)= accel(startInd:startInd+399,i);
    vel(:,i)= -cumtrapz(acc(:,i));
%     kk=(vel(1,i)-vel(400,i))/400;
%     for j=1:400
%         vel(j,i)=vel(j,i)+j*kk;
%     end;
    vel(:,i)=vel(:,i)*29/max(vel(:,i));
    [maxVel(i) maxInd(i)] = max(vel(:,i));
end;
meanVel=sum(vel,2)/length(total_trials);

[minVI,minI]= min(maxInd);
[maxVI,maxI]= max(maxInd);
[meanV,meanI]=max(meanVel);

figure(3)
t=(1:400)*5;
plot(t,meanVel);
hold on;
plot(t,vel(:,minI), 'g-');
hold on;
plot(t,vel(:,maxI), 'r-');
hold on;
plot([5*meanI,5*meanI],[-3,meanV],'b-');
xlabel('Time(ms)');

fprintf('minV  minI   maxV  maxI   meanV meanI\n');
fprintf('%6.2f %6.2f  %6.2f %6.2f %6.2f %6.2f\n', maxVel(minI),minVI,maxVel(maxI),maxVI,meanV,meanI);

minDelay=minVI*5-1000
maxDelay=maxVI*5-1000
meanDelay=meanI*5-1000

minSyn=min(firstSyn)
maxSyn=max(firstSyn)

figure(4)
subplot(3,1,1),hist(visualBegintime);
subplot(3,1,2),hist(firstSyn-visualBegintime);
subplot(3,1,3),hist(maxInd*5-1000);


   
% velocity_hist = 0;
% velocity_hist_foraft = 0;
% % for j = 1 : 1000
% %     velocity_hist = velocity_hist + ( data.eye_data(5,j,:)-offset_y ) * 0.2;
% %  %   velocity_hist = velocity_hist + ( data.eye_data(5,j,:)-mean(data.eye_data(5,1:1000,:)) ) * 0.2;
% %     velocity_int(j,:) = velocity_hist;   
% %     velocity_hist_foraft = velocity_hist_foraft + ( data.eye_data(6,j,:)-offset_y_foraft ) * 0.2;
% %     velocity_int_foraft(j,:) = velocity_hist_foraft;  
% % end
% for j = 1 : 200        
%     velocity_hist = velocity_hist + ( data.eye_data(5,300+23+j,:)-mean(data.eye_data(5,201+23:250+23,:)) ) * 0.2; %use first 250ms
%   %  velocity_hist = velocity_hist + (
%   %  data.eye_data(5,200+23+j,:)-mean(data.eye_data(5,201+23:600+23,:)) ) *0.2; % use 2s
%   %velocity_hist = velocity_hist + ( data.eye_data(5,300+23+j,:)-mean(data.eye_data(5,301+23:500+23,:)) ) * 0.2; % use 1s
%     velocity_int(j,:) = velocity_hist;  
% end
% 
% for run = 1:1
%     
%     for i = 1 : length(spike_rates)
%    %     resp_accelerometer(1,i) = sum(velocity_int(1+(run-1)*4:20+(run-1)*4,i+BegTrial-1))/200;   % normalize velocity
%         resp_accelerometer(1,i) = sum(velocity_int(1:200,i+BegTrial-1))/200;
%    %   resp_accelerometer(1,i) = sum(velocity_int(101:300,i+BegTrial-1))/200;
%   %      resp_accelerometer_foraft(1,i) = sum(velocity_int_foraft(300:500,i+BegTrial-1))/200;    % normalize velocity 
%    %     accelerometer_raw(:,i) = squeeze(data.eye_data(5,300:500,i))-mean(squeeze(data.eye_data(5,300:500,i)));
%     end
% 
%     % select = find( heading == unique_heading(5) & (stim_type == 1)  ) ;
%     % for i = 1:length(spike_rates)
%     %
%     %%%%%%%%%%%%%%%%%%%%%%%%Edited to diagnose centrifuge moog. (AB-2/05/07)
%     % % for i = 1:length(stim_type) % when using heading discrim 2I only!
%     % for i = 1:length(spike_rates)
%     %     a6(i,:) = data.eye_data(6,:, i);  
%     % end
%     % accel(1,:) = mean(a6(:,:));
%     % figure
%     % plot(accel(1,:))
%     % title('acceleration trace');
%     % % dlmwrite('accel_155.txt',accel(1,:))
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     % % consider syncpulse
%     % select = 
%     % temp = find(data.spike_data(2,:,i)==1);
%     % temp2 = temp(temp>1000 ); 
%     % syncpulse_index(i) = temp2(1); % align with first syncpulse   
%     % for i=length(spike_rates)
%     %     if syncpulse_index(i)<1200 %sometimes it is strange that syncpulse does not happen in a proper time
%     %         a6(i,:) = data.eye_data(6,round(syncpulse_index(i)/5) : round(syncpulse_index(i)/5)+399,i); % spike data has 5 times higher resolution than eye data
%     %     else
%     %        a6(i,:) = data.eye_data(6,201 : 600, i);  
%     %     end
%     % end
%     % accel(1,:) = median(a6(select,:));
% 
%     %determine for each trial whether monkey chooses leftward(target1) or rightward(tarket2)
%     LEFT = 1;
%     RIGHT = 2;
%     for i=1 : length(spike_rates)
%         temp = data.event_data(1,:,i+BegTrial-1);
%         events = temp(temp>0);  % all non-zero entries
%         if (sum(events == IN_T1_WIN_CD) > 0)
%             choice(i) = RIGHT;
%         elseif (sum(events == IN_T2_WIN_CD) > 0)
%             choice(i) = LEFT;
%         else
%             disp('Neither T1 or T2 chosen.  This should not happen!.  File must be bogus.');
%         end
%     end
% 
%     % psychometric dataset
%     psycho_correct = [];
%     fit_data_psycho = [];
%     N_obs = [];
%     for i = 1:length(unique_heading)
%         for k = 1:length(unique_stim_type)
%              trials_p =logical( (heading == unique_heading(i)) & (stim_type == unique_stim_type(k)) ) ;
%              % make 'S' curve by using the rightward choice for y-axis
%              correct_trials = (trials_p & (choice == RIGHT) );
%              psycho_correct(k,i) = 1*sum(correct_trials) / sum(trials_p); 
%              fit_data_psycho_cum{k}(i, 1) = unique_heading( i );  
%              fit_data_psycho_cum{k}(i, 2) = psycho_correct(k,i);
%              fit_data_psycho_cum{k}(i, 3) = sum(trials_p); 
%          end
%     end
% 
%     resp_heading = [];
%     Z_Spikes = resp_accelerometer;
%     % z-score data for later cp analysis across headings
%     for k = 1:length(unique_stim_type)
%         for i = 1:length(unique_heading)
%             select =logical( (heading == unique_heading(i)) & (stim_type == unique_stim_type(k)) ) ;  
%             z_dist = resp_accelerometer(select);
%             z_dist = (z_dist - mean(z_dist))/std(z_dist);
%             Z_Spikes(select) = z_dist;
%         end
%     end
% 
%     % now group neuronal data into two groups according to monkey's choice
%     for k = 1:length(unique_stim_type)    % notice that the condition is double than disc_heading  
%         select0 =logical( (heading == 0) & (stim_type == 1) ) ;
%         selectleft = logical( (heading == unique_heading(1)) & (stim_type == 1) ) ;
%         selectright = logical( (heading == unique_heading(end)) & (stim_type == 1) ) ;
%         for i = 1:length(unique_heading)
%             select =logical( (heading == unique_heading(i)) & (stim_type == unique_stim_type(k)) ) ;   
%             selectves =logical( (heading == unique_heading(i)) & (stim_type == unique_stim_type(1)) ) ;        
%             resp{k,i} = resp_accelerometer(select); 
%             % calculate firing rate of each trial
%     %         for j = 1 : repetition; 
%     %             spike_temp = resp_accelerometer(select);   
%     %             resp_heading{k}(i, j) = spike_temp( j );           
%     %         end
%             resp_mat{k}(i) = mean(resp{k,i});  % the mean firing rate for each heading 
%             resp_mat_err{k}(i) = std(resp{k,i}) / sqrt(repetition);
%             % calculate CP, group data based on monkey's choice 
%             resp_left_choose{k,i} = resp_accelerometer(select & (choice == LEFT) );
%             resp_right_choose{k,i} = resp_accelerometer(select & (choice == RIGHT) );
%             if (length(resp_left_choose{k,i}) <= 3) | (length(resp_right_choose{k,i}) <= 3)   % make sure each condition has at least 3 data values
%           %  if (length(resp_left_choose{k,i}) / length(resp{k,i}) <0.25) |  (length(resp_left_choose{k,i}) / length(resp{k,i}) >0.75)  
%                 Z_Spikes(select) = 9999;   % similar to NaN, just make a mark
%             end 
%             if  isempty( resp_accelerometer(selectves & (choice == LEFT)) ) ==0 
%                 resp_hist_leftchoice(i,1:121) = hist( resp_accelerometer(selectves & (choice == LEFT) ), -3:0.05:3);
%             else
%                 resp_hist_leftchoice(i,1:121) = 0;
%             end
%             if isempty( resp_accelerometer(selectves & (choice == RIGHT)) ) ==0 
%                resp_hist_rightchoice(i,1:121) = hist( resp_accelerometer(selectves & (choice == RIGHT) ),  -3:0.05:3);
%             else
%                resp_hist_rightchoice(i,1:121) = 0;
%             end
% %             accelerometer_raw_mean(i,k) = mean(mean(accelerometer_raw(:,select)));
% %             accelerometer_raw_std(i,k) = std(std(accelerometer_raw(:,select)));
%         end  
% 
%         % now across all data 
%         resp_left_all{k} = Z_Spikes(  (stim_type == unique_stim_type(k)) & (choice == LEFT) & (Z_Spikes~=9999) ); 
%         resp_right_all{k} = Z_Spikes(  (stim_type == unique_stim_type(k)) & (choice == RIGHT) & (Z_Spikes~=9999) ); 
%         resp_left{k} = Z_Spikes(  (stim_type == unique_stim_type(k)) & (heading == 0) & (choice == LEFT) & (Z_Spikes~=9999) ); 
%         resp_right{k} = Z_Spikes(  (stim_type == unique_stim_type(k)) & (heading == 0) & (choice == RIGHT) & (Z_Spikes~=9999) ); 
%         resp_all{k} = Z_Spikes(  (stim_type == unique_stim_type(k)) & (Z_Spikes~=9999) ); 
%     end
% 
%     %-------------------------------------------------------------------------
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
%                %  Neuro_correct{k}(i) =  rocN( resp{k,(1+length(unique_heading))/2},resp{k,i},100 ); % compare to the 0 heading condition, which is straght ahead
%                  Neuro_correct{k}(i) =  rocN( resp{k,length(unique_heading)-i+1},resp{k,i},100 );    % antineuron model
%              else
%                %  Neuro_correct{k}(i) =  rocN( resp{k,(1+length(unique_heading))/2}, resp{k,(i+1)},100 ); % compare to the 0 heading condition, which is straght ahead
%                  Neuro_correct{k}(i) =  rocN( resp{k,length(unique_heading)-i}, resp{k,(i+1)},100 );   % antineuron model  
%              end
%              if  resp_mat{k}(1) < resp_mat{k}(end)
%                  Neuro_correct{k}(i) = 1 - Neuro_correct{k}(i);             
%              end  
%          end
%     end
%     % choice probability
%     for k = 1 : length(unique_stim_type)
%         for i = 1 : length(unique_heading)  
%             if (length(resp_left_all{k}) > 3) & (length(resp_right_all{k}) > 3)
%             %if  (length(resp_left_all{k}) / length(resp_all{k}) >0.25) |  (length(resp_left_all{k}) / length(resp_all{k}) <0.75)
%                 CP_all{k}(run) = rocN( resp_left_all{k},resp_right_all{k},100 );
%             else
%                 CP_all{k}(run) = NaN;
%             end
%             if (length(resp_left_choose{k,i}) > 3) & (length(resp_right_choose{k,i}) > 3)
%                CP{k}(run,i) = rocN( resp_left_choose{k,i},resp_right_choose{k,i},100 );
%             else
%                 CP{k}(run,i) = NaN;
%             end
%             if  resp_mat{k}(1) < resp_mat{k}(end)
%                 CP_all{k}(run) = 1 - CP_all{k}(run);
%                 CP{k}(run,i) = 1 - CP{k}(run,i);
%             end
%         end
%     end    
% 
% % % %--------------------------------------------------------------------------
% % %do permutation to test the significance of CP_all{k}, re-calculate CP 2000 times
% % perm_num = 1000;
% % Z_Spikes_perm = Z_Spikes;
% % bin = 0.005;
% % x_bin = 0 : bin : 1;
% % for k = 1:length(unique_stim_type)    
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
% %         resp_left_perm{k} = Z_Spikes_perm(  (stim_type == unique_stim_type(k)) & (heading == 0) & (choice == LEFT) & (Z_Spikes~=9999) ); 
% %         resp_right_perm{k} = Z_Spikes_perm(  (stim_type == unique_stim_type(k)) & (heading == 0) & (choice == RIGHT) & (Z_Spikes~=9999) ); 
% %         
% %         if  (length(resp_left_all{k}) > 3) & (length(resp_right_all{k}) > 3) 
% %             CP_all_perm{k}(n) = rocN( resp_left_all_perm{k}, resp_right_all_perm{k},100 );
% %         else
% %             CP_all_perm{k}(n) = NaN; 
% %         end
% %         if  (length(resp_left{k}) > 3) & (length(resp_right{k}) > 3) 
% %             CP_perm{k}(n) = rocN( resp_left_perm{k}, resp_right_perm{k},100 );
% %         else
% %             CP_perm{k}(n) = NaN; 
% %         end
% %         if  resp_mat{k}(1) < resp_mat{k}(end)  
% %             CP_all_perm{k}(n) = 1 - CP_all_perm{k}(n);     
% %             CP_perm{k}(n) = 1 - CP_perm{k}(n); 
% %         end  
% %     end
% %     % now calculate p value or significant test
% %     if (length(resp_left_all{k}) > 3) & (length(resp_right_all{k}) > 3) 
% %         hist_perm(k,:) = hist( CP_all_perm{k}(:), x_bin );  % for permutation
% %         bin_sum = 0;
% %         n = 0;
% %         while ( n < (CP_all{k}/bin) )
% %              n = n+1;
% %              bin_sum = bin_sum + hist_perm(k, n);
% %              if CP_all{k} > 0.5                  % note it's two tail test
% %                 p{k} = 2*(perm_num - bin_sum)/ perm_num;    % calculate p value for CP_all
% %              else
% %                 p{k} = 2* bin_sum / perm_num;
% %              end
% %         end
% %     else
% %         p{k} = NaN;
% %     end 
% %     p_0{k} = length(find(CP{k}>CP_perm{k}))/1000;
% % end
% % 
%     %%%%%% use Wichman's MLE method to estimate threshold and bias
%     for k = 1:length(unique_stim_type)        
% %         wichman_psy = pfit(fit_data_psycho_cum{k},'plot_opt','no plot','shape','cumulative gaussian','n_intervals',1,'FIX_LAMBDA',0.001,'sens',0,'compute_stats','false','verbose','false');  
% %         Thresh_psy{k} = wichman_psy.params.est(2);
% %         Bias_psy{k} = wichman_psy.params.est(1);
% %         psy_perf{k} = [wichman_psy.params.est(1),wichman_psy.params.est(2)];
%         fit_data_neuro_cum{k}(:,1) = unique_heading(unique_heading~=0);
%         fit_data_neuro_cum{k}(:,2) = Neuro_correct{k}(:);
% %         wichman_neu = pfit(fit_data_neuro_cum{k},'plot_opt','no plot','shape','cumulative gaussian','n_intervals',1,'FIX_LAMBDA',0.001,'sens',0,'compute_stats','false','verbose','false');
% %         Thresh_neu{k} = wichman_neu.params.est(2);
% %         Bias_neu{k} = wichman_neu.params.est(1);
% %         neu_perf{k} = [wichman_neu.params.est(1),wichman_neu.params.est(2)];
%         
%         [bb, tt] = cum_gaussfit_max1(fit_data_neuro_cum{k});
%         if tt>300 | tt<0
%             tt=300;
%         end
%         threshold_accelerometer{k}(run) = tt;
%     end
% 
% end
% 
% % %plot psychometric and neurometric function here
% % h{1} = 'bo';  f{1} = 'b-';  g{1} = 'bo-';
% % h{2} = 'rd';  f{2} = 'r-';  g{2} = 'rd-';
% % h{3} = 'gs';  f{3} = 'g-';  g{3} = 'gs-';
% % figure(10);
% % set(10,'Position', [5,25, 980,650], 'Name', 'psycho_neurometic function');
% % orient landscape;
% % %plot psychometric function
% % axes('position',[0.05 0.3, 0.26 0.45]);
% % title('psychometric');
% % legend_txt = [];
% % for k = 1:length(unique_stim_type)
% %     xi = min(unique_heading) : 0.1 : max(unique_heading);   
% %     beta = [0, 1.0];
% %  %   plot data in logarithmical space instead of linspace
% %     plot(unique_heading, psycho_correct(k,:), h{k}, xi, cum_gaussfit(psy_perf{k}, xi),  f{k} );
% %     xlabel('Heading Angles');   
% %     ylim([0,1]);
% %     ylabel('Rightward Choices');
% %     hold on;
% %     legend_txt{k*2-1} = [num2str(unique_stim_type(k))];
% %     legend_txt{k*2} = [''];
% % %    also fit data with weibull function
% % %    [psycho_alpha(k) psycho_beta(k)]= weibull_fit(fit_data_psycho{k});
% % end
% % %plot neurometric function
% % axes('position',[0.36 0.3, 0.26 0.45]);
% % title('neuroometric');
% % for k = 1:length(unique_stim_type)
% %     neu_heading = unique_heading(unique_heading~=0);
% %     xi = min(unique_heading) : 0.1 : max(unique_heading); 
% %     plot(neu_heading, Neuro_correct{k}(:), h{k},  xi, cum_gaussfit(neu_perf{k}, xi),  f{k} );
% %     xlabel('Heading Angles');   
% %     ylim([0,1]);
% %     hold on;
% % %    betafit_ne_cutt(k)=betafit_ne_cut{k}(2);
% % %    also fit data with weibull function
% % %     [neuro_alpha(k) neuro_beta(k)]= weibull_fit(fit_data_neuro{k});
% % %     [neuro_alpha_cut(k) neuro_beta_cut(k)]= weibull_fit(fit_data_neuro_cut{k});  % tick the most outside datapoint out
% % end
% % 
% % %%%%%%  neurological raw data based on firing rate instead of ROC
% % axes('position',[0.7 0.3, 0.26 0.45]);
% % title('firing rate');
% % for k = 1:length(unique_stim_type)
% %     errorbar(unique_heading, resp_mat{k}(:), resp_mat_err{k}(:),g{k} );
% %     xlabel('Heading Angle (deg)');
% %     ylabel('Firing rate(spikes/s)');   
% %     xlim([min(unique_heading),max(unique_heading)]);
% %     hold on;
% % end
% % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % output some text of basic parameters in the figure
% % axes('position',[0.05,0.8, 0.9,0.15] );
% % xlim( [0,100] );
% % ylim( [2,10] );
% % text(0, 10, FILE);
% % text(20,10,'coherence =');
% % text(30,10,num2str(unique_motion_coherence) );
% % text(45,10,'repetition =');
% % text(55,10,num2str(repetition) ); 
% % text(5,8, 'Psy: u      threshold        err           Neu:u         threshold         err              CP        p');
% % text(0,8, 'stim');
% % for k = 1:length(unique_stim_type)
% %     text(0,8-k, num2str(unique_stim_type(k)));
% %     text(5,8-k,num2str(Bias_psy{k} ));
% %     text(12,8-k,num2str(Thresh_psy{k} ));
% % %    text(20,8-k,num2str(psy_boot(k)) );
% %     text(30,8-k,num2str(Bias_neu{k} ));
% %     text(40,8-k,num2str(Thresh_neu{k} ));
% %   %  text(50,8-k,num2str(neu_boot(k)) );
% %     text(50,8-k,num2str(CP_all{k}') ); 
% %     text(60,8-k,num2str(p{k}') );   
% %  %   text(53,8-k,num2str(CP{k}(3:end-2)) );  % always show the middle angles, not the 24 and 8 
% % end
% % axis off;
% % %--------------------------------------------------------------------------
% %Also, write out some summary data to a cumulative summary file
% sprint_txt = ['%s']; 
% for i = 1 : 5000 % this should be large enough to cover all the data that need to be exported
%      sprint_txt = [sprint_txt, ' %4.3f'];    
% end
% buff1 = sprintf( sprint_txt, FILE, threshold_accelerometer{1} );
% outfile1 = ['Z:\Users\Yong\accelerometerThreshold_250ms.dat'];
% printflag = 0;
% if (exist(outfile1, 'file') == 0)   % file does not yet exist
%     printflag = 1;
% end
% fid = fopen(outfile1, 'a');
% if (printflag)
%     fprintf(fid, 'FILE\t');
%     fprintf(fid, '\r\n');
% end
% fprintf(fid, '%s', buff1);
% fprintf(fid, '\r\n');
% fclose(fid);
% 
% buff2 = sprintf(sprint_txt, FILE, CP_all{1}, CP{1}(1,5) );
% outfile2 = ['Z:\Users\Yong\accelerometerCP_250ms.dat'];
% printflag = 0;
% if (exist(outfile2, 'file') == 0)   % file does not yet exist
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

%---------------------------------------------------------------------------------------
return;