function FiringRate(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

Path_Defs;
ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP

%get the column of values for azimuth and elevation and stim_type
temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG);
temp_azimuth = data.moog_params(AZIMUTH,:,MOOG);
temp_elevation = data.moog_params(ELEVATION,:,MOOG);

%now, get the firing rates for all the trials 
temp_spike_rates = data.spike_rates(SpikeChan, :);  
temp_spike_data = data.spike_data(1,:);   % spike rasters
%temp_spike_rates = FIR_Filter(temp_spike_rates, 20, 100, 'high', 20, 0);

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
spike_rates = temp_spike_rates(~null_trials & select_trials);

unique_azimuth = munique(azimuth');
unique_elevation = munique(elevation');
unique_stim_type = munique(stim_type');

condition = stim_type;
h_title{1}='Vestibular';
h_title{2}='Visual';
h_title{3}='Combined';
unique_condition = munique(condition');

% use spike_data to compute mean firing rate
Discard_trials = find(null_trials==1 | trials <BegTrial | trials >EndTrial);
for i = 1 : length(Discard_trials)
    temp_spike_data( 1, ((Discard_trials(i)-1)*5000+1) :  Discard_trials(i)*5000 ) = 99;
end
spike_data(1,:) = temp_spike_data( 1, find(temp_spike_data(1,:)~=99) );
spike_data(1, find(spike_data>10) ) = 1; % something is absolutely wrong 

% calculate spontaneous firing rate
spon_found = find(null_trials==1); 
spon_resp = mean(temp_spike_rates(spon_found));
%repetition = floor( length(spike_rates) / (26*length(unique_stim_type)) ); % take minimum repetition
pc=0;
for k=1: length(unique_stim_type)
    for j=1:length(unique_elevation)
        for i=1:length(unique_azimuth)        
            select = find( (azimuth==unique_azimuth(i)) & (elevation==unique_elevation(j)) & (stim_type==unique_stim_type(k)) );
            if (sum(select) > 0) 
                pc=pc+1;
                trialrepeat(pc) = length(select);
            end
        end
    end
end
repetition = min(trialrepeat)

StartEventBin(1)=996;
windowlength = 50; % 50 ms
for w = 1:1
windowlength = windowlength + (w-1)*50;
%count = 0;
for s = 1 : 1
    % replace spike_rates with spike_data based on analize window set    
%     for ss =  1 : length(spike_rates) % ss marks the index of trial
% %         spike_rates(ss) = sum( spike_data(1,StartEventBin(1)+windowlength*(s-1)+5000*(ss-1) : StartEventBin(1)+windowlength+windowlength*(s-1)+5000*(ss-1)) ) ; % 996~3006 every 200ms
%           spike_rates(ss) = sum( spike_data(1,StartEventBin(1)+500+5000*(ss-1) : StartEventBin(1)+500+1000+5000*(ss-1)) );%Enlarge the window
%     end    

    %% ADD CODE HERE FOR PLOTTING
    if length(unique_azimuth)==10
        unique_azimuth0=[0:45:315]';
    else
        unique_azimuth0=unique_azimuth;
    end
    resp_mat = [];    
    for k=1: length(unique_condition)        
 %       for j=1:length(unique_elevation) 
        resp_trial_temp =[];
        resp_trial_group = [];
        count = 0;
        %for j=1:length(unique_elevation) 
        for j=1:1 
            for i=1:length(unique_azimuth)
         %       select = logical( (azimuth==unique_azimuth(i)) & (elevation==unique_elevation(j)) & (condition==unique_condition(k)) );
                select = logical( (azimuth==unique_azimuth(i)) & (elevation==0) & (condition==unique_condition(k)) );%AC 12-11-2007
                if (sum(select) > 0)  
                    count=count+1;
                    spike_temp = spike_rates(select);
                    resp_mat_anova_horizontal{k}(1:repetition,i) = spike_temp(1:repetition);
                    resp_trial_temp = [resp_trial_temp, spike_temp];
                    resp_trial_group_temp =[];
                    resp_trial_group_temp(1:length(spike_temp)) = i;
                    resp_trial_group = [resp_trial_group,resp_trial_group_temp]; 
                    
                    resp_mat(k, j, i) = mean(spike_rates(select));
                    resp_std(k, j, i) = std(spike_rates(select));
                    resp_std_root(k, j, i) = std(sqrt(spike_rates(select)));
                    resp_mat_count(k,count) = resp_mat(k, j, i);
                    resp_mat_std_count(k,count) = resp_std(k, j, i);
                    
                    % z-score data                  
                    z_dist = spike_rates(select);
                    if std(z_dist)~=0 % there are cases that all values are 0 for a certain condition, e.g. m2c73r1, visual condition
                       z_dist = (z_dist - mean(z_dist))/std(z_dist);
                    else
                        z_dist = 0;
                    end
                    Z_Spikes(select) = z_dist;  
                    
%                 else
%                     resp_mat(k, j, i) = resp_mat(k,j,1);
%                     resp_std(k, j, i) = resp_std(k,j,1);
%                     resp_std_root(k, j, i) = resp_std_root(k, j, 1);
                end
            end        
        end
        resp_trial{k}(:, 1) = resp_trial_temp;
        resp_trial{k}(:, 2) = resp_trial_group;
        P_anova_horizontal(k) = anovan(resp_trial{k}(:,1),{resp_trial{k}(:,2)},'display','off');  
        [p_anova, table, stats] = anova1(resp_mat_anova_horizontal{k},[],'off');
        P_anova_horizontal2(k) = p_anova;
    end

%     for k = 1 : length(unique_condition) 
%         resp_mat_26(k,:) = [resp_mat(k,1,1),squeeze(resp_mat(k,2,:))',squeeze(resp_mat(k,3,:))',squeeze(resp_mat(k,4,:))',resp_mat(k,5,1)];
%         resp_mat_std_26(k,:) = [resp_std(k,1,1),squeeze(resp_std(k,2,:))',squeeze(resp_std(k,3,:))',squeeze(resp_std(k,4,:))',resp_std(k,5,1)];
%     end

    % calculate anova
    % trials_per_rep = (length(unique_azimuth)*length(unique_elevation)-2*(length(unique_azimuth)-1)) * length(unique_condition) + 1;
    % repetitions = floor( (EndTrial-(BegTrial-1)) / trials_per_rep);
    %repetition = floor( (EndTrial-(BegTrial-1)) / 79);

    % resp_mat_anova = [];
    % for k=1: length(unique_condition) 
    %     n=0;
    %     for i=1:length(unique_azimuth)
    % %         for j=1:length(unique_elevation)
    % %             select_rep = find( azimuth==unique_azimuth(i) & elevation==unique_elevation(j) & condition==unique_condition(k) );
    % %             select_rep_horizontal = find( azimuth==unique_azimuth(i) & elevation==0 & condition==unique_condition(k) );
    % %             if (length(select_rep) > 0)    
    % %                 n = n+1;            
    % %                 for q=1:repetitions
    % %                    resp_mat_anova{k}(q,n) = spike_rates(select_rep(q));
    % %                    resp_mat_anova_horizontal{k}(q,n) = spike_rates(select_rep_horizontal(q));
    % %                 end
    % %             end
    % %         end
    %         trial_select = logical( (azimuth==unique_azimuth(i)) & elevation==0 & (condition==unique_condition(k)) );
    %         for jj = 1 : repetitions; 
    %             spike_temp = spike_rates(trial_select);   
    %             resp_horizontal_trial{k}(jj, i) = spike_temp( jj );  
    %         end        
    %    end
    % %    [p_anova, table, stats] = anova1(resp_mat_anova{k},[],'off');
    % %    P_anova(k) = p_anova;
    % %    [p_anova, table, stats] = anova1(resp_mat_anova_horizontal{k},[],'off');
    % %    P_anova_horizontal(k) = p_anova;
    % end
    % deal different conditions
    % if length(unique_condition) ==1
    %    visfind = 1;
    % elseif length(unique_condition) ==2
    %    visfind = 2; 
    % elseif length(unique_condition) ==3
    %    visfind = 2; 
    % end
    %ves(1:26) = [resp_mat(1,1,1),squeeze(resp_mat(1,2,:))',squeeze(resp_mat(1,3,:))',squeeze(resp_mat(1,4,:))',resp_mat(1,5,1)];
    % vis(1:26) = [resp_mat(visfind,1,1),squeeze(resp_mat(visfind,2,:))',squeeze(resp_mat(visfind,3,:))',squeeze(resp_mat(visfind,4,:))',resp_mat(visfind,5,1)];
    % com(1:26) = [resp_mat(3,1,1),squeeze(resp_mat(3,2,:))',squeeze(resp_mat(3,3,:))',squeeze(resp_mat(3,4,:))',resp_mat(3,5,1)];
    %ves_std(1:26) = [resp_std(1,1,1),squeeze(resp_std(1,2,:))',squeeze(resp_std(1,3,:))',squeeze(resp_std(1,4,:))',resp_std(1,5,1)];
    % vis_std(1:26) = [resp_std(visfind,1,1),squeeze(resp_std(visfind,2,:))',squeeze(resp_std(visfind,3,:))',squeeze(resp_std(visfind,4,:))',resp_std(visfind,5,1)];
    %com_std(1:26) = [resp_std(3,1,1),squeeze(resp_std(3,2,:))',squeeze(resp_std(3,3,:))',squeeze(resp_std(3,4,:))',resp_std(3,5,1)];
    %vesout = resp_mat(1, 3, :);
    %visout = resp_mat(2, 3, :);
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %Also, write out some summary data to a cumulative summary file
%     sprint_txt = ['%s']; 
%     for i = 1 : 500 % this should be large enough to cover all the data that need to be exported
%          sprint_txt = [sprint_txt, ' %4.3f'];    
%     end
%     % buff = sprintf(sprint_txt, FILE, Z_Spikes );
%     %buff = sprintf(sprint_txt, FILE, spon_resp, vis, vis_std, P_anova(visfind)  );
%     %buff = sprintf(sprint_txt, FILE,  resp_mat(1,3,:),resp_mat(2,3,:),resp_std(1,3,:),resp_std(2,3,:),resp_std_root(1,3,:),resp_std_root(2,3,:),P_anova_horizontal(1),P_anova_horizontal(2)  );
%     buff = sprintf(sprint_txt, FILE, resp_mat_26(:,:), resp_mat_std_26(:,:).^2 );
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
%     %outfile = ['Z:\Data\Tempo\Batch Files\Tunde\dat_files\Yong_3D_horizontalplane.dat'];
%     %outfile = ['C:\Aihua\z_TempOutputs\VesOnly_Aihua.dat'];
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
%---------------------------------------------------------------------------------------
eye_x_left_temp_temp(:,:) = data.eye_data(1,:,~null_trials & select_trials );
eye_y_left_temp_temp(:,:) = data.eye_data(2,:,~null_trials & select_trials );
eye_x_right_temp_temp(:,:) = data.eye_data(3,:,~null_trials & select_trials );
eye_y_right_temp_temp(:,:) = data.eye_data(4,:,~null_trials & select_trials );

eye_x_left_temp(:,:) = eye_x_left_temp_temp(:, elevation==0);
eye_y_left_temp(:,:) = eye_y_left_temp_temp(:, elevation==0);
eye_x_right_temp(:,:) = eye_x_right_temp_temp(:, elevation==0);
eye_y_right_temp(:,:) = eye_y_right_temp_temp(:, elevation==0);
dim1 = size(eye_x_left_temp);
for i=1:dim1(2)
    eyeleftmaxmin(i) = abs( max(eye_x_left_temp(1:600,i))-min(eye_x_left_temp(1:600,i)) );
    eyerightmaxmin(i) = abs( max(eye_x_right_temp(1:600,i))-min(eye_x_right_temp(1:600,i)) );
end
for i=1:200
    eye_x_left(i,:) = eye_x_left_temp(i+322,:)-mean(eye_x_left_temp(222:322,:));
    eye_y_left(i,:) = eye_y_left_temp(i+322,:)-mean(eye_y_left_temp(222:322,:));
    eye_x_right(i,:) = eye_x_right_temp(i+322,:)-mean(eye_x_right_temp(222:322,:));
    eye_y_right(i,:) = eye_y_right_temp(i+322,:)-mean(eye_y_right_temp(222:322,:));
    deviation_left(i,:) = sqrt(eye_x_left(i,:).^2+eye_y_left(i,:).^2);
    deviation_right(i,:) = sqrt(eye_x_right(i,:).^2+eye_y_right(i,:).^2);
end
if median(eyeleftmaxmin)<median(eyerightmaxmin) % left eye has no signal
    deviation(1,:) = median(deviation_right(:,:),2);
    eye_pos = deviation_right(:,:);
else
    deviation(1,:) = median(deviation_left(:,:),2);
    eye_pos = deviation_left(:,:);
end
dim = size(eye_pos);
count = 0;
for j=1:dim(2)    
    gauss2 = normpdf(1:1:200,100,1); % smooth data at a SD of 5ms (1 point)
    eye_temp = conv(eye_pos(:,j), gauss2);
    eye_pos_smooth(:,j) = eye_temp(100:end-ceil(200/2)); 
    eye_vel(:,j) = abs(diff(eye_pos_smooth(:,j)))*200; 
    saccadefind = find( eye_vel(:,j)>10 );
    if length(saccadefind)>=1 
       tempfind = find(diff(saccadefind)~=1);
       if length(tempfind)>=1
           saccade_count = length(tempfind)+1; % the number of microscaddes            
           tempfind2 = [0 tempfind' length(saccadefind)];
           for n=1:saccade_count
               peak_vel(n+count) = max( eye_vel(saccadefind(tempfind2(n)+1):saccadefind(tempfind2(n+1)),j) );
           end
       else
           saccade_count = 1;           
           peak_vel(1+count) = max(eye_vel(:,j));
       end
    else
        saccade_count = 0;        
    end   
    count = count+saccade_count; 
end
saccade_rate=count/dim(2);
saccade_vel = hist(log10(peak_vel), 1:0.05:2)/count;

%Also, write out some summary data to a cumulative summary file
sprint_txt = ['%s']; 
for i = 1 : 1000 % this should be large enough to cover all the data that need to be exported
     sprint_txt = [sprint_txt, ' %4.3f'];    
end
% buff = sprintf(sprint_txt, FILE, length(unique_stim_type), rr );
% outfile = ['Z:\Users\Yong\variancetomeancorrelation50ms.dat'];
% buff = sprintf(sprint_txt, FILE, resp_mat_count, resp_mat_std_count.^2 );

%buff = sprintf(sprint_txt,FILE,8,repetition,resp_mat_anova_horizontal{1}(1:repetition,:),resp_mat_anova_horizontal{2}(1:repetition,:) );
%buff = sprintf(sprint_txt,FILE,resp_mat_count(1,:),resp_mat_count(2,:),resp_mat_count(3,:));
%buff = sprintf(sprint_txt,resp_mat_count(1,:),resp_mat_count(2,:),resp_mat_count(3,:));
buff = sprintf(sprint_txt, FILE, saccade_rate, deviation, saccade_vel);
outfile = ['Z:\Users\Yong\eyetrain.dat'];

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
end

return;