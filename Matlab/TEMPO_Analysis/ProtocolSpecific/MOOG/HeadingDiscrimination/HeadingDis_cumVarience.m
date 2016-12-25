%-----------------------------------------------------------------------------------------------------------------------
%-- psychometric and neurometric function for heading discrimination task
%--	07/16/04 GY
%-----------------------------------------------------------------------------------------------------------------------

function HeadingDis_cumVarience(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin,StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

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
%temp_spike_rates = FIR_Filter(temp_spike_rates, 20, 100, 'high', 20, 0); %filter
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this part deals with for whatever reason, some trials are absolutely
% wrong, like >1000 spikes/s, makes no sense. tick this bad trial out,
% replace with one of the other data within same condition, this is similar
% to bootstrap. But in order to keep every running the same, use the mean
% of the rest instead. this happen only very rarely!!!
bad_tri = find(temp_spike_rates > 1000 );
if length(bad_tri) > 0
    for i=1:length(bad_tri)
		same_tri = find( (stim_type == stim_type(bad_tri(i))) & (heading == heading(bad_tri(i))) ); % find other trials within same condition
		rest_tri = same_tri(same_tri ~= bad_tri(i));
		spike_rates(bad_tri(i)) = mean( spike_rates(rest_tri) );
    end
end

% use spike_data to compute mean firing rate
Discard_trials = find( trials <BegTrial | trials >EndTrial);
for i = 1 : length(Discard_trials)
    temp_spike_data( 1, ((Discard_trials(i)-1)*5000+1) :  Discard_trials(i)*5000 ) = 9999;
end
spike_data(1,:) = temp_spike_data( 1, find(temp_spike_data(1,:)~=9999) );
spike_data(1, find(spike_data>10) ) = 1; % something is absolutely wrong 

one_repetition = length(unique_heading)*length(unique_stim_type);
repetition = floor( length(spike_rates)/one_repetition ); % take minimum repetition

StartEventBin(1)=996;
outlier =[];
windowlength = 50; % 50 ms
for w=1:1
windowlength = windowlength + (w-1)*50;
%count = 0;
for s = 1 : 1
    % replace spike_rates with spike_data based on analize window set    
%     for ss =  1 : length(spike_rates) % ss marks the index of trial
% %         spike_rates(ss) = sum( spike_data(1,StartEventBin(1)+windowlength*(s-1)+5000*(ss-1) : StartEventBin(1)+windowlength+windowlength*(s-1)+5000*(ss-1)) ) ; % 996~3006 every 200ms
%        if w==1 
%           spike_rates(ss) = sum( spike_data(1,StartEventBin(1)+500+5000*(ss-1) : StartEventBin(1)+500+1000+5000*(ss-1)) );%Enlarge the window
%        else
%           spike_rates(ss) = sum( spike_data(1,StartEventBin(1)+750+5000*(ss-1) : StartEventBin(1)+750+500+5000*(ss-1)) );
%        end
%     end      
    %repetition = repetition -1;
    resp_heading = [];
    % for cell m2c384r2, only 44 repeats, one trial is missing for k=2,i=8;
    %repetition = 44;
    for c = 1:length(unique_motion_coherence)
        count = 0;
        for k = 1:length(unique_stim_type)    % notice that the condition is double than disc_heading    
            for i = 1:length(unique_heading)
                count = count+1;
                select =logical( (heading == unique_heading(i)) & (stim_type == unique_stim_type(k)) & (motion_coherence == unique_motion_coherence(c))) ;  
                act_found = find( select==1 );
                resp{k,i} = spike_rates(select);   
        %         % calculate firing rate of each trial  
        %         for j = 1 : repetition; 
        %             spike_temp = spike_rates(select);   
        %             resp_trial{k}(j, i) = spike_temp(j);           
        %         end  
                resp_mat{c}(k,i) = mean(resp{k,i});  % the mean firing rate for each heading 
                resp_mat_std{c}(k,i) = std(resp{k,i}); % std    
                resp_mat_count(w,count) = resp_mat{c}(k,i);
                resp_mat_std_count(w,count) = resp_mat_std{c}(k,i);

                % z-score data for spike count correlation analysis
                z_dist = spike_rates(select);
                if std(z_dist)~=0 % there are cases that all values are 0 for a certain condition, e.g. m2c73r1, visual condition
                   z_dist = (z_dist - mean(z_dist))/std(z_dist);
                else
                    z_dist = 0;
                end
                Z_Spikes(select) = z_dist;   

            end

        %     pp = polyfit([1,2,3,4,5], resp_mat{k}(3:7), 1); % fitting a line to the 5 data pts about straight ahead
        %     slope(k) = pp(1);
        %     offset(k) = pp(2);
        %     ave_std(k) = mean( sqrt(resp_mat_varience{k}(3:7)) ); % the average std for the five data pts I'm interested in
        %     
        %     [RRR,PPP] = corrcoef(resp_mat{k}(:), resp_mat_varience{k}(:));
        %     r(k) = RRR(1,2);
        %     p(k) = PPP(1,2);

        %     pp = polyfit(resp_mat{k}(3:7), resp_mat_varience{k}(:), 1);
        %     slope(k) = pp(1);
        %     offset(k) = pp(2);
        %     [RRR,PPP] = corrcoef(resp_mat{k}(:), resp_mat_varience{k}(:));
        %     r(k) = RRR(1,2);
        %     p(k) = PPP(1,2);
        %   pp(k) = anova1(resp_trial{k},'','off');
        end
    end
%     z_spikes_window(s, :) = Z_Spikes;
%     % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %Also, write out some summary data to a cumulative summary file
%     sprint_txt = ['%s']; 
%     for i = 1 : 500 % this should be large enough to cover all the data that need to be exported
%          sprint_txt = [sprint_txt, ' %4.3f'];    
%     end
%     % buff = sprintf( sprint_txt, FILE, repetition, resp_trial{1}(1,:),resp_trial{1}(2,:),resp_trial{1}(3,:),resp_trial{1}(4,:),resp_trial{1}(5,:),resp_trial{1}(6,:),resp_trial{1}(7,:),resp_trial{1}(8,:),resp_trial{1}(9,:), ... 
%     %                                                resp_trial{2}(1,:),resp_trial{2}(2,:),resp_trial{2}(3,:),resp_trial{2}(4,:),resp_trial{2}(5,:),resp_trial{2}(6,:),resp_trial{2}(7,:),resp_trial{2}(8,:),resp_trial{2}(9,:), ...
%     %                                                resp_trial{3}(1,:),resp_trial{3}(2,:),resp_trial{3}(3,:),resp_trial{3}(4,:),resp_trial{3}(5,:),resp_trial{3}(6,:),resp_trial{3}(7,:),resp_trial{3}(8,:),resp_trial{3}(9,:) );
%     %buff = sprintf(sprint_txt, FILE, slope, ave_std, resp_mat{1}, resp_mat_std{1}, resp_mat{2}, resp_mat_std{2},resp_mat{3}, resp_mat_std{3} );
%     buff = sprintf(sprint_txt, FILE, resp_mat(:,:), resp_mat_std(:,:).^2 );
%     % buff = sprintf(sprint_txt, FILE, repetition,resp_trial{1}(1, :),resp_trial{1}(2, :),resp_trial{1}(3, :),resp_trial{1}(4, :),resp_trial{1}(5, :),resp_trial{1}(6, :),resp_trial{1}(7, :),resp_trial{1}(8, :),resp_trial{1}(9, :), ...
%     %                                             resp_trial{2}(1, :),resp_trial{2}(2, :),resp_trial{2}(3, :),resp_trial{2}(4, :),resp_trial{2}(5, :),resp_trial{2}(6, :),resp_trial{2}(7, :),resp_trial{2}(8, :),resp_trial{2}(9, :), ...
%     %                                             resp_trial{3}(1, :),resp_trial{3}(2, :),resp_trial{3}(3, :),resp_trial{3}(4, :),resp_trial{3}(5, :),resp_trial{3}(6, :),resp_trial{3}(7, :),resp_trial{3}(8, :),resp_trial{3}(9, :) );
% 
%     if s==1
%        outfile = ['Z:\Users\Yong\variancetomean1.dat'];
%     elseif s==2
%         outfile = ['Z:\Users\Yong\variancetomean2.dat'];
%     elseif s==3
%         outfile = ['Z:\Users\Yong\variancetomean3.dat'];
%     elseif s==4
%         outfile = ['Z:\Users\Yong\variancetomean4.dat'];
%     elseif s==5
%         outfile = ['Z:\Users\Yong\variancetomean5.dat'];
%     elseif s==6
%         outfile = ['Z:\Users\Yong\variancetomean6.dat'];
%     elseif s==7
%         outfile = ['Z:\Users\Yong\variancetomean7.dat'];
%     elseif s==8
%         outfile = ['Z:\Users\Yong\variancetomean8.dat'];
%     else
%         outfile = ['Z:\Users\Yong\variancetomean9.dat'];
%     end
%     
%     printflag = 0;
%     if (exist(outfile, 'file') == 0)   % file does not yet exist
%         printflag = 1;
%     end
%     fid = fopen(outfile, 'a');
%     if (printflag)
%         fprintf(fid, 'FILE\t');
%         fprintf(fid, '\r\n');
%     end
%     fprintf(fid, '%s', buff);
%     fprintf(fid, '\r\n');
%     fclose(fid);  
    z_spikes_window(s,:)=Z_Spikes;    
end


% %z-score across for serial temporal correlation
% for k = 1:length(unique_stim_type) % separated among stimuli conditions
%     select = logical( stim_type == unique_stim_type(k) );
%     z_spikes_window_stim = [];    
%     z_spikes_window_stim(:,:) = z_spikes_window(:,select);
%     count = 0;   
%     count_stim = 0;
%     for i = 1:2000/windowlength 
%         for j = 1:2000/windowlength             
%            count = count+1;
%            count_stim = count_stim+1;
%            index = find(z_spikes_window(i,:)<=3 & z_spikes_window(i,:)>=-3 & z_spikes_window(j,:)<=3 & z_spikes_window(j,:)>=-3);
%            index_stim = find(z_spikes_window_stim(i,:)<=3 & z_spikes_window_stim(i,:)>=-3 & z_spikes_window_stim(j,:)<=3 & z_spikes_window_stim(j,:)>=-3);
%            if sum(z_spikes_window(i,index))~=0 & sum(z_spikes_window(j,index))~=0
%               [r p] = corrcoef( z_spikes_window(i,index), z_spikes_window(j,index) );
%               rr(count) = r(1,2);   
%            else
%               rr(count)=0;
%            end
%            if sum(z_spikes_window_stim(i,index_stim))~=0 & sum(z_spikes_window_stim(j,index_stim))~=0
%                [r p] = corrcoef( z_spikes_window_stim(i,index_stim), z_spikes_window_stim(j,index_stim) );
%                rr_stim(k,count_stim) = r(1,2); 
%            else
%                rr_stim(k,count_stim) = 0;
%            end
%         end
%     end
% end
end
%Also, write out some summary data to a cumulative summary file
sprint_txt = ['%s']; 
for i = 1 : 1000 % this should be large enough to cover all the data that need to be exported
     sprint_txt = [sprint_txt, ' %4.3f'];    
end
% buff = sprintf(sprint_txt, FILE, length(unique_stim_type), rr );
% outfile = ['Z:\Users\Yong\variancetomeancorrelation50ms.dat'];
% buff = sprintf(sprint_txt, FILE, resp_mat_count, resp_mat_std_count.^2 );
% if w==1
%     buff = sprintf(sprint_txt, FILE, resp_mat_count(1,:), resp_mat_std_count(1,:).^2 );
%    outfile = ['Z:\Users\Yong\variancetomean_middle1000ms.dat'];
% else
%     buff = sprintf(sprint_txt, FILE, resp_mat_count(2,:), resp_mat_std_count(2,:).^2 );
%    outfile = ['Z:\Users\Yong\variancetomean_middle500ms.dat'];
% end
buff = sprintf(sprint_txt, FILE, resp_mat{1}(1,:), resp_mat_std{1}(1,:), resp_mat{1}(2,:), resp_mat_std{1}(2,:), resp_mat{1}(3,:), resp_mat_std{1}(3,:));
outfile = ['Z:\Users\Yong\fisherMeanVarianceAllProtocol000ms.dat'];

printflag = 0;
if (exist(outfile, 'file') == 0)   % file does not yet exist
    printflag = 1;
end
fid = fopen(outfile, 'a');
if (printflag)
    fprintf(fid, 'FILE\t');
    fprintf(fid, '\r\n');
end
fprintf(fid, '%s', buff);
fprintf(fid, '\r\n');
fclose(fid);  

return;