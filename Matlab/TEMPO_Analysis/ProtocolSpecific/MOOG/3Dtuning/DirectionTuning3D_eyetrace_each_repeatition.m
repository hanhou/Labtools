function DirectionTuning3D_eyetrace(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

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
temp_spike_data = data.spike_data(1,:);   % spike rasters

null_trials = logical( (temp_azimuth == data.one_time_params(NULL_VALUE)) );
%now, remove trials from direction and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:length(temp_azimuth);
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );
azimuth = temp_azimuth(~null_trials & select_trials);
elevation = temp_elevation(~null_trials & select_trials);
stim_type = temp_stim_type(~null_trials & select_trials);
amplitude = temp_amplitude(~null_trials & select_trials);
spike_rates = temp_spike_rates(~null_trials & select_trials);

unique_azimuth = munique(azimuth');
unique_elevation = munique(elevation');
unique_stim_type = munique(stim_type');
unique_amplitude = munique(amplitude');

repeat = floor( length(temp_spike_rates) / (length(unique_stim_type)*(length(unique_azimuth)*(length(unique_elevation)-2)+2)+1) );

% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % calculate mean
% % for k = 1 : length(unique_stim_type)
% %     select_exp = find( (((temp_azimuth==90)&(temp_elevation==0)) | ((temp_azimuth==45)&(temp_elevation==0)) | ((temp_azimuth==135)&(temp_elevation ...
% %                         ==0)) | ((temp_azimuth==90)&(temp_elevation==45)) | ((temp_azimuth==90)&(temp_elevation==-45)))  & (temp_stim_type==unique_stim_type(k)) );
% %     select_con = find( (((temp_azimuth==270)&(temp_elevation==0)) | ((temp_azimuth==225)&(temp_elevation==0)) | ((temp_azimuth==315)&(temp_elevation ...
% %                         ==0)) | ((temp_azimuth==270)&(temp_elevation==45)) | ((temp_azimuth==270)&(temp_elevation==-45)))  & (temp_stim_type==unique_stim_type(k)) );
% %     group1=[];
% %     group2=[];
% %     for i = 1 : length(select_exp)
% %         verg_exp{k}(1,i) = mean( data.eye_data(1,301:500,select_exp(i)) - data.eye_data(3,301:500,select_exp(i)) );
% %         group1 = [group1,1]; % for anovan
% %     end
% %     for j = 1 : length(select_con)
% %         verg_con{k}(1,j) = mean( data.eye_data(1,301:500,select_con(j)) - data.eye_data(3,301:500,select_con(j)) );
% %         group2 = [group2,2]; % for anovan
% %     end
% %     group3 = [group1,group2];
% %     group = {group3};
% %     p(k) = anovan([verg_exp{k}(1,:),verg_con{k}(1,:)], group,'linear',3,'off');
% %     close;
% % end
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % find expansion and contraction trials
% for k = 1 : length(unique_stim_type)
%      select_exp = find( (((temp_azimuth==90)&(temp_elevation==0)) | ((temp_azimuth==45)&(temp_elevation==0)) | ((temp_azimuth==135)&(temp_elevation ...
%                          ==0)) | ((temp_azimuth==90)&(temp_elevation==45)) | ((temp_azimuth==90)&(temp_elevation==-45)))  & (temp_stim_type==unique_stim_type(k)) );
%      select_con = find( (((temp_azimuth==270)&(temp_elevation==0)) | ((temp_azimuth==225)&(temp_elevation==0)) | ((temp_azimuth==315)&(temp_elevation ...
%                          ==0)) | ((temp_azimuth==270)&(temp_elevation==45)) | ((temp_azimuth==270)&(temp_elevation==-45)))  & (temp_stim_type==unique_stim_type(k)) );
%      % only use straight forward and backward point
%      select_f = find( (temp_azimuth==90)&(temp_elevation==0)&(temp_stim_type==unique_stim_type(k)) );
%      select_b = find( (temp_azimuth==270)&(temp_elevation==0)&(temp_stim_type==unique_stim_type(k)) );     
% %     for i = 1 : 400
% %         verg_exp_mean(k,i) = mean( (data.eye_data(1,i+200,select_exp)-data.eye_data(3,i+200,select_exp)) );
% %       %  verg_exp_std(k,i) = std( data.eye_data(1,i+200,select_exp)-data.eye_data(3,i+200,select_exp) );
% %         verg_con_mean(k,i) = mean( (data.eye_data(1,i+200,select_con)-data.eye_data(3,i+200,select_con)) );
% %       %  verg_con_std(k,i) = std( data.eye_data(1,i+200,select_con)-data.eye_data(3,i+200,select_con) );       
% %     end
% %     verg_exp(k,:) = verg_exp_mean(k,:)-mean(verg_exp_mean(k,:));
% %     verg_con(k,:) = verg_con_mean(k,:)-mean(verg_con_mean(k,:));
%     for t = 1 : repeat*5
%         temp_exp = select_exp(t);
%         temp_con = select_con(t);
%         resp_exp_trail(k,t) = temp_spike_rates(temp_exp);
%         resp_con_trail(k,t) = temp_spike_rates(temp_con);
%         % take the mean again based on the middle 1 sec
%         verg_exp_trail(k,t) = mean( (data.eye_data(1,301:500,select_exp(t))-data.eye_data(3,301:500,select_exp(t))) ); 
%         verg_con_trail(k,t) = mean( (data.eye_data(1,301:500,select_con(t))-data.eye_data(3,301:500,select_con(t))) ); 
%     end 
%     for t1 = 1 : repeat
%         temp_f = select_f(t1);
%         temp_b = select_b(t1);
%         resp_f_trail(k,t1) = temp_spike_rates(temp_f);
%         resp_b_trail(k,t1) = temp_spike_rates(temp_b);
%         % take the mean again based on the middle 1 sec
%         verg_f_trail(k,t1) = mean( (data.eye_data(1,301:500,select_f(t1))-data.eye_data(3,301:500,select_f(t1))) ); 
%         verg_b_trail(k,t1) = mean( (data.eye_data(1,301:500,select_b(t1))-data.eye_data(3,301:500,select_b(t1))) );         
%     end
%     % group data in all azimuth within the horizontal plane
%     for aa = 1 : length(unique_azimuth) % 8 azimuth
%         select_azimuth = find( (temp_elevation==0) & (temp_stim_type==unique_stim_type(k)) & (temp_azimuth==unique_azimuth(aa)) );
%         for tt = 1 : repeat
%             verg_azi(k,aa,tt) = mean( (data.eye_data(1,301:500,select_azimuth(tt))-data.eye_data(3,301:500,select_azimuth(tt))) );
%             verg_azi_left(k,aa,tt) = mean( data.eye_data(1,301:500,select_azimuth(tt)) ) - mean( data.eye_data(1,211:300,select_azimuth(tt)) );
%             verg_azi_right(k,aa,tt) = mean( data.eye_data(3,301:500,select_azimuth(tt)) ) - mean( data.eye_data(3,211:300,select_azimuth(tt)) );
% %            verg_azi(k,aa,tt) = verg_azi(k,aa,tt) - mean(mean ( data.eye_data(1,201:300,:) - data.eye_data(3,201:300,:) )); % subtract DC offset
%         end        
%     end
% end
% azimuth_verg(:,1:8) = squeeze(verg_azi(1,:,:))';
% azimuth_verg(:,9:16) = squeeze(verg_azi(2,:,:))';
% azimuth_verg(:,17:24) = squeeze(verg_azi(3,:,:))';
% azimuth_verg_LR(:,1:8) = squeeze(verg_azi_left(1,:,:))';
% azimuth_verg_LR(:,9:16) = squeeze(verg_azi_left(2,:,:))';
% azimuth_verg_LR(:,17:24) = squeeze(verg_azi_left(3,:,:))';
% azimuth_verg_LR(:,25:32) = squeeze(verg_azi_right(1,:,:))';
% azimuth_verg_LR(:,33:40) = squeeze(verg_azi_right(2,:,:))';
% azimuth_verg_LR(:,41:48) = squeeze(verg_azi_right(3,:,:))';
% pl{1}='b-'; pl{2}='r-'; pl{3}='g-'; pl{4}='k-'; pl{5}='y-'; 
% pr{1}='bo'; pr{2}='ro'; pr{3}='go'; pr{4}='ko'; pr{5}='yo'; 
% select_in = find( (temp_elevation==0) );
% figure(2);
% set(2,'Position', [5,15 980,650], 'Name', '3D Direction Tuning');
% x=1:400;
% for i=1:repeat*8*3
%     plot(data.eye_data(1,201:600,select_in(i))', x', 'b-');   
%     set(gca, 'ydir' , 'reverse');    
%     hold on;
%     plot(data.eye_data(3,201:600,select_in(i))', x', 'r-');
%     set(gca, 'ydir' , 'reverse');
% end
% figure(3);
% set(3,'Position', [15,5 980,650], 'Name', '3D Direction Tuning');
% x=1:400;
% for i=1:repeat*8*3
%     plot(data.eye_data(1,201:600,select_in(i))'-data.eye_data(3,201:600,select_in(i))', x', 'b-');   
%     set(gca, 'ydir' , 'reverse');    
%     hold on;
% %     plot(data.eye_data(3,201:600,i)'-mean(data.eye_data(3,201:300,i)), x', 'r-');
% %     set(gca, 'ydir' , 'reverse');
% end
% for i = 1: repeat*8*3
%     aa(i,1) = mean(data.eye_data(1,201:300,select_in(i)) - data.eye_data(3,201:300,select_in(i)) );
%     vv(i,1) = mean(data.eye_data(1,301:500,select_in(i)) - data.eye_data(3,301:500,select_in(i)) );
% end

%--------------------------------------------------------------------------
% for Katzu's analysis
for k = 1 : length(unique_stim_type)
    for i = 1 : length(unique_azimuth)
        for j = 1 : length(unique_elevation)
            select = find( temp_azimuth==unique_azimuth(i) & temp_elevation==unique_elevation(j) & temp_stim_type==unique_stim_type(k) );
            if sum(select)>0
                for jj = 1 : repeat % for convenience with sacrefice of some of the trials
               %For Que and Zebulon
%                     offset_x = mean( data.eye_data(3,201:300,select(jj)) ); % horizontal
%                     offset_y = mean( data.eye_data(4,201:300,select(jj)) ); % vertical
%                     resp_x{k,i,j}(jj,:) = data.eye_data(3,201:600,select(jj)) - offset_x;  % horizontal   
%                     resp_y{k,i,j}(jj,:) = data.eye_data(4,201:600,select(jj)) - offset_y;  % horizontal   
                %For Azrael
                    offset_x = mean( data.eye_data(1,201:300,select(jj)) ); % horizontal
                    offset_y = mean( data.eye_data(2,201:300,select(jj)) ); % vertical
                    resp_x{k,i,j}(jj,:) = data.eye_data(1,201:600,select(jj)) - offset_x;  % horizontal   
                    resp_y{k,i,j}(jj,:) = data.eye_data(2,201:600,select(jj)) - offset_y;  % horizontal   
                end
            else
                resp_x{k,i,j}(:,:) = resp_x{k,1,j}(:,:);
                resp_y{k,i,j}(:,:) = resp_y{k,1,j}(:,:);
            end
        end
    end
    resp_x_up{k}(:,:) = resp_x{k,1,1}(:,:);     resp_y_up{k}(:,:) = resp_y{k,1,1}(:,:);
    resp_x_down{k}(:,:) = resp_x{k,1,5}(:,:);   resp_y_down{k}(:,:) = resp_y{k,1,5}(:,:);
    resp_x_left{k}(:,:) = resp_x{k,5,3}(:,:);   resp_y_left{k}(:,:) = resp_y{k,5,3}(:,:);
    resp_x_right{k}(:,:) = resp_x{k,1,3}(:,:);  resp_y_right{k}(:,:) = resp_y{k,1,3}(:,:);  
    % reshape the matrix for output text file
    re_x_up{k} = reshape(resp_x_up{k}(:,:)',1,repeat*400);
    re_y_up{k} = reshape(resp_y_up{k}(:,:)',1,repeat*400);
    re_x_down{k} = reshape(resp_x_down{k}(:,:)',1,repeat*400);
    re_y_down{k} = reshape(resp_y_down{k}(:,:)',1,repeat*400);
    re_x_left{k} = reshape(resp_x_left{k}(:,:)',1,repeat*400);
    re_y_left{k} = reshape(resp_y_left{k}(:,:)',1,repeat*400);
    re_x_right{k} = reshape(resp_x_right{k}(:,:)',1,repeat*400);
    re_y_right{k} = reshape(resp_y_right{k}(:,:)',1,repeat*400);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% output to text file
sprint_txt = ['%s'];
for i = 1 : 400*repeat*10
    sprint_txt = [sprint_txt, ' %1.2f'];    
end
% if you want to save 'vestibular', select 1, 'visual', select 2

% buff= sprintf(sprint_txt, FILE, repeat, re_x_up{1}(:,:), re_y_up{1}(:,:), re_x_down{1}(:,:), re_y_down{1}(:,:), ...
%                           re_x_left{1}(:,:), re_y_left{1}(:,:),re_x_right{1}(:,:),re_y_right{1}(:,:) );
buff= sprintf(sprint_txt, FILE, repeat, re_x_up{2}(:,:), re_y_up{2}(:,:), re_x_down{2}(:,:), re_y_down{2}(:,:), ...
                          re_x_left{2}(:,:), re_y_left{2}(:,:),re_x_right{2}(:,:),re_y_right{2}(:,:) );


% outfile = [BASE_PATH 'ProtocolSpecific\MOOG\3Dtuning\Eye_tra_Zebulon_ves.dat'];
% outfile = [BASE_PATH 'ProtocolSpecific\MOOG\3Dtuning\Eye_tra_Zebulon_vis.dat'];

% outfile = [BASE_PATH 'ProtocolSpecific\MOOG\3Dtuning\Eye_tra_Azrael_ves.dat'];
outfile = [BASE_PATH 'ProtocolSpecific\MOOG\3Dtuning\Eye_tra_Azrael_vis.dat'];


% outfile = [BASE_PATH 'ProtocolSpecific\MOOG\3Dtuning\Eye_tra_Que_ves.dat'];
% outfile = [BASE_PATH 'ProtocolSpecific\MOOG\3Dtuning\Eye_tra_Que_vis.dat'];

printflag = 0;
if (exist(outfile, 'file') == 0)    %file does not yet exist
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

