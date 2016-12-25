% function DirectionTuning3D_eyetrace(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
% 
% TEMPO_Defs;
% Path_Defs;
% ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % for figure writing , if running batch, comment
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % From here RVOR/pursuit
% % if size(data.eye_data,1)>6
% %     LEFT_EYE_1_2=9;
% %     RIGHT_EYE_3_4=10;
% % else
% %     LEFT_EYE_1_2=7;
% %     RIGHT_EYE_3_4=8;
% % end
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %************************************************************************%
% % %plot Vertical vs. Horizontal
% % switch (data.eye_flag)
% %     case (LEFT_EYE_1_2)
% %         Eye_Select='Left Eye'
% %         Hor=1;        Ver=2;
% %     case(RIGHT_EYE_3_4)
% %         Eye_Select='Right Eye'
% %         Hor=3;        Ver=4;
% % end
% % % 
% 
% temp_azimuth = data.moog_params(AZIMUTH,:,MOOG);
% temp_elevation = data.moog_params(ELEVATION,:,MOOG);
% temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG);
% temp_spike_rates = data.spike_rates(SpikeChan, :); 
% temp_total_trials = data.misc_params(OUTCOME, :);
% temp_spike_data = data.spike_data(1,:);   % spike rasters
% % temp_fp_rotate = data.moog_params(FP_ROTATE,:,MOOG);
% 
% % null_trials = logical( (temp_azimuth == data.one_time_params(NULL_VALUE)) );
% null_trials = logical( (temp_elevation == data.one_time_params(NULL_VALUE)) );
% %now, remove trials from direction and spike_rates that do not fall between BegTrial and EndTrial
% % trials = 1:length(temp_azimuth);
% trials = 1:length(temp_elevation);
% select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );
% azimuth = temp_azimuth(~null_trials & select_trials);
% elevation = temp_elevation(~null_trials & select_trials);
% stim_type = temp_stim_type(~null_trials & select_trials);
% spike_rates = temp_spike_rates(~null_trials & select_trials);
% % fp_rotate = temp_fp_rotate(~null_trials & select_trials);
% 
% unique_azimuth = munique(azimuth');
% unique_elevation = munique(elevation');
% unique_stim_type = munique(stim_type');
% % unique_fp_rotate = munique(fp_rotate');
% 
% % repeat = floor( length(temp_spike_rates) / (length(unique_stim_type)*(length(unique_fp_rotate)*(length(unique_elevation)-2)+2)+1) );
% %Yong's bad Aihua's good for calculate repeat
% 
% trials_per_rep = (length(unique_azimuth)*length(unique_elevation)-14) * length(unique_stim_type) + 1;
% repeat = floor( (EndTrial-(BegTrial-1)) / trials_per_rep);
% 
% for k = 1 : length(unique_stim_type)
%     for i = 1 : length(unique_azimuth)
%         for j = 1 : length(unique_elevation)
%             select = find( temp_azimuth==unique_azimuth(i) & temp_elevation==unique_elevation(j) & temp_stim_type==unique_stim_type(k) );
%             select_up = find( temp_azimuth==0 & temp_elevation==-90 & temp_stim_type==1 );
%             if sum(select)>0
%                 for jj = 1 : repeat % for convenience with sacrefice of some of the trials
% %                Right Eye Que_Rotation, old_Que_c1_c71(Tra), Zeblon,
% %                Lothar, Que_c290_c333(Tra), Que_Labyrinthectomy,
% %                Zebulon_Laby
%                     offset_x = mean( data.eye_data(3,201:300,select(jj)) ); % horizontal
%                     offset_y = mean( data.eye_data(4,201:300,select(jj)) ); % vertical
%                     resp_x{k,i,j}(jj,:) = data.eye_data(3,201:600,select(jj)) - offset_x;  % horizontal   
%                     resp_y{k,i,j}(jj,:) = data.eye_data(4,201:600,select(jj)) - offset_y;  % horizontal 
% %                Left Eye Azrael, Lothar(after..),
% %                new_Que_c86__butCh.3,4 are O.K., too
% %                     offset_x = mean( data.eye_data(1,201:300,select(jj)) ); % horizontal
% %                     offset_y = mean( data.eye_data(2,201:300,select(jj)) ); % vertical
% %                     resp_x{k,i,j}(jj,:) = data.eye_data(1,201:600,select(jj)) - offset_x;  % horizontal   
% %                     resp_y{k,i,j}(jj,:) = data.eye_data(2,201:600,select(jj)) - offset_y;  % horizontal
% %                Eye select (for figure wrinting)
% %                     offset_x = mean( data.eye_data(Hor,201:300,select(jj)) ); % horizontal
% %                     offset_y = mean( data.eye_data(Ver,201:300,select(jj)) ); % vertical
% %                     resp_x{k,i,j}(jj,:) = data.eye_data(Hor,201:600,select(jj)) - offset_x;  % horizontal   
% %                     resp_y{k,i,j}(jj,:) = data.eye_data(Ver,201:600,select(jj)) - offset_y;  % vertical
%                 end
%             else
%                 resp_x{k,i,j}(:,:) = resp_x{k,1,j}(:,:);%??? +-90 elevation
%                 resp_y{k,i,j}(:,:) = resp_y{k,1,j}(:,:);
%             end
%         end
%     end
%     
%     % resp_x{k,i,1}(:,:) i=1 means stim_type = vestibular, if you want to
%     % visual turn to resp_x{k,2,1}(:,:)
%     
%     resp_x_up{k}(:,:) = resp_x{k,1,1}(:,:);     resp_y_up{k}(:,:) = resp_y{k,1,1}(:,:);
%     resp_x_down{k}(:,:) = resp_x{k,1,5}(:,:);   resp_y_down{k}(:,:) = resp_y{k,1,5}(:,:);
%     resp_x_left{k}(:,:) = resp_x{k,5,3}(:,:);   resp_y_left{k}(:,:) = resp_y{k,5,3}(:,:);
%     resp_x_right{k}(:,:) = resp_x{k,1,3}(:,:);  resp_y_right{k}(:,:) = resp_y{k,1,3}(:,:);  
% %     resp_x_up{k}(:,:) = resp_x{k,2,5}(:,:);     resp_y_up{k}(:,:) = resp_y{k,2,5}(:,:);
% %     resp_x_down{k}(:,:) = resp_x{k,2,1}(:,:);   resp_y_down{k}(:,:) = resp_y{k,2,1}(:,:);
% %     resp_x_left{k}(:,:) = resp_x{k,2,3}(:,:);   resp_y_left{k}(:,:) = resp_y{k,2,3}(:,:);
% %     resp_x_right{k}(:,:) = resp_x{k,2,7}(:,:);  resp_y_right{k}(:,:) = resp_y{k,2,7}(:,:);  
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% %%%%%   Yong's Eye analysis matlab file  No more use
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % reshape the matrix for output text file
%     %% 400points-400-400-400...400*5repeats=2000points only for fixation
%     
%     
% %     re_x_up{k} = reshape(resp_x_up{k}(:,:)',1,repeat*400);
% %     re_y_up{k} = reshape(resp_y_up{k}(:,:)',1,repeat*400);
% %     re_x_down{k} = reshape(resp_x_down{k}(:,:)',1,repeat*400);
% %     re_y_down{k} = reshape(resp_y_down{k}(:,:)',1,repeat*400);
% %     re_x_left{k} = reshape(resp_x_left{k}(:,:)',1,repeat*400);
% %     re_y_left{k} = reshape(resp_y_left{k}(:,:)',1,repeat*400);
% %     re_x_right{k} = reshape(resp_x_right{k}(:,:)',1,repeat*400);
% %     re_y_right{k} = reshape(resp_y_right{k}(:,:)',1,repeat*400);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %     mean for 5 repeatition for figures
%     re_x_up{k} = mean(resp_x_up{k}(:,:));
%     re_y_up{k} = mean(resp_y_up{k}(:,:));
%     re_x_down{k} = mean(resp_x_down{k}(:,:));
%     re_y_down{k} = mean(resp_y_down{k}(:,:));
%     re_x_left{k} = mean(resp_x_left{k}(:,:));
%     re_y_left{k} = mean(resp_y_left{k}(:,:));
%     re_x_right{k} = mean(resp_x_right{k}(:,:));
%     re_y_right{k} = mean(resp_y_right{k}(:,:));
% %     
% 
% end
% % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% %For figure writing, in roder to batch comment!!!!
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% 
% 
%    title1 = 'up';
%     title2 = 'down';
%     title3 = 'left';
%     title4 = 'right';
%     
%     
%     figure(2)
% 
% subplot(4,1,1)
% plot(re_x_up{1},'r.');
%     hold on;
%     xlim( [1, 400] );
%     set(gca, 'XTickMode','manual');
%     set(gca, 'xtick',[1,100,200,300,400]);
%     set(gca, 'xticklabel','0|0.5|1|1.5|2'); 
% %     ylim([-20, 20]);
%     ylabel('(deg)');
%     title(['Eye Position /  ',title1]);
% 
%     plot(re_y_up{1},'b.');
%     hold off;
% 
% 
% subplot(4,1,2)
% plot(re_x_down{1},'r.');
%     hold on;
%     xlim( [1, 400] );
%     set(gca, 'XTickMode','manual');
%     set(gca, 'xtick',[1,100,200,300,400]);
%     set(gca, 'xticklabel','0|0.5|1|1.5|2'); 
% %     ylim([-20, 20]);
%     ylabel('(deg)');
%     title(['Eye Position /  ',title2]);
% 
%     plot(re_y_down{1},'b.');
%     hold off;
% 
% subplot(4,1,3)
% plot(re_x_left{1},'r.');
%     hold on;
%     xlim( [1, 400] );
%     set(gca, 'XTickMode','manual');
%     set(gca, 'xtick',[1,100,200,300,400]);
%     set(gca, 'xticklabel','0|0.5|1|1.5|2'); 
% %     ylim([-20, 20]);
%     ylabel('(deg)');
%     title(['Eye Position /  ',title3]);
% 
%     plot(re_y_left{1},'b.');
%     hold off;
%     
% subplot(4,1,4)
% plot(re_x_right{1},'r.');
%     hold on;
%     xlim( [1, 400] );
%     set(gca, 'XTickMode','manual');
%     set(gca, 'xtick',[1,100,200,300,400]);
%     set(gca, 'xticklabel','0|0.5|1|1.5|2'); 
% %     ylim([-20, 20]);
%     ylabel('(deg)');
%     title(['Eye Position /  ',title4]);
% 
%     plot(re_y_right{1},'b.');
%     hold off;
%     
%  figure(3)
% 
% subplot(4,1,1)
% plot(re_x_up{2},'r.');
%     hold on;
%     xlim( [1, 400] );
%     set(gca, 'XTickMode','manual');
%     set(gca, 'xtick',[1,100,200,300,400]);
%     set(gca, 'xticklabel','0|0.5|1|1.5|2'); 
% %     ylim([-20, 20]);
%     ylabel('(deg)');
%     title(['Eye Position /  ',title1]);
% 
%     plot(re_y_up{2},'b.');
%     hold off;
% 
% 
% subplot(4,1,2)
% plot(re_x_down{2},'r.');
%     hold on;
%     xlim( [1, 400] );
%     set(gca, 'XTickMode','manual');
%     set(gca, 'xtick',[1,100,200,300,400]);
%     set(gca, 'xticklabel','0|0.5|1|1.5|2'); 
% %     ylim([-20, 20]);
%     ylabel('(deg)');
%     title(['Eye Position /  ',title2]);
% 
%     plot(re_y_down{2},'b.');
%     hold off;
% 
% subplot(4,1,3)
% plot(re_x_left{2},'r.');
%     hold on;
%     xlim( [1, 400] );
%     set(gca, 'XTickMode','manual');
%     set(gca, 'xtick',[1,100,200,300,400]);
%     set(gca, 'xticklabel','0|0.5|1|1.5|2'); 
% %     ylim([-20, 20]);
%     ylabel('(deg)');
%     title(['Eye Position /  ',title3]);
% 
%     plot(re_y_left{2},'b.');
%     hold off;
%     
% subplot(4,1,4)
% plot(re_x_right{2},'r.');
%     hold on;
%     xlim( [1, 400] );
%     set(gca, 'XTickMode','manual');
%     set(gca, 'xtick',[1,100,200,300,400]);
%     set(gca, 'xticklabel','0|0.5|1|1.5|2'); 
% %     ylim([-20, 20]);
%     ylabel('(deg)');
%     title(['Eye Position /  ',title4]);
% 
%     plot(re_y_right{2},'b.');
%     hold off;
%     
%     
%     
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%                         output to text file
% sprint_txt = ['%s'];
% for i = 1 : 400*repeat*10
%     sprint_txt = [sprint_txt, ' %1.2f'];    
% end
% 
% 
% % if you want to save 'world', select 2, 'pursuit', select 3
% 
% 
% buff= sprintf(sprint_txt, FILE, repeat, re_x_up{1}(:,:), re_y_up{1}(:,:), re_x_down{1}(:,:), re_y_down{1}(:,:), ...
%                           re_x_left{1}(:,:), re_y_left{1}(:,:),re_x_right{1}(:,:),re_y_right{1}(:,:) );
% % buff= sprintf(sprint_txt, FILE, repeat, re_x_up{2}(:,:), re_y_up{2}(:,:), re_x_down{2}(:,:), re_y_down{2}(:,:), ...
% %                           re_x_left{2}(:,:), re_y_left{2}(:,:),re_x_right{2}(:,:),re_y_right{2}(:,:) );
% % buff= sprintf(sprint_txt, FILE, repeat, re_x_up{3}(:,:), re_y_up{3}(:,:), re_x_down{3}(:,:), re_y_down{3}(:,:), ...
% %                           re_x_left{3}(:,:), re_y_left{3}(:,:),re_x_right{3}(:,:),re_y_right{3}(:,:) );
% %                       
%                       
% % outfile = [BASE_PATH 'ProtocolSpecific\MOOG\RVOR_Pursuit\Que_Pursuit_World.dat'];
% % outfile = [BASE_PATH 'ProtocolSpecific\MOOG\3Dtuning\AzraelEye_Tra_Vet.dat'];
% % outfile = [BASE_PATH 'ProtocolSpecific\MOOG\3Dtuning\AzraelEye_Tra_Vis.dat'];
% 
% % outfile = [BASE_PATH 'ProtocolSpecific\MOOG\3Dtuning\QueEye_Tra_Vet.dat'];
% % outfile = [BASE_PATH 'ProtocolSpecific\MOOG\3Dtuning\ZebulonEye_Tra_Vet.dat'];
% 
% % outfile = [BASE_PATH 'ProtocolSpecific\MOOG\3Dtuning\QueLaby_Tra_Vet.dat'];
% outfile = [BASE_PATH 'ProtocolSpecific\MOOG\3Dtuning\ZebulonLaby_Tra_Vet.dat'];
% 
% printflag = 0;
% if (exist(outfile, 'file') == 0)    %file does not yet exist
%     printflag = 1;
% end
% fid = fopen(outfile, 'a');
% if (printflag)
%     fprintf(fid, 'FILE\t');
%     fprintf(fid, '\r\n');
% end
% fprintf(fid, '%s', buff);
% fprintf(fid, '\r\n');
% fclose(fid);
% 
% 
% 
% %%%%%                         output to text file
% sprint_txt2 = ['%s'];
% for i = 1 : 400*repeat*10
%     sprint_txt2 = [sprint_txt2, ' %1.2f'];    
% end
% 
% 
% % if you want to save 'world', select 2, 'pursuit', select 3
% 
% 
% % buff= sprintf(sprint_txt2, FILE, repeat, re_x_up{1}(:,:), re_y_up{1}(:,:), re_x_down{1}(:,:), re_y_down{1}(:,:), ...
% %                           re_x_left{1}(:,:), re_y_left{1}(:,:),re_x_right{1}(:,:),re_y_right{1}(:,:) );
% buff= sprintf(sprint_txt, FILE, repeat, re_x_up{2}(:,:), re_y_up{2}(:,:), re_x_down{2}(:,:), re_y_down{2}(:,:), ...
%                           re_x_left{2}(:,:), re_y_left{2}(:,:),re_x_right{2}(:,:),re_y_right{2}(:,:) );
% % buff= sprintf(sprint_txt, FILE, repeat, re_x_up{3}(:,:), re_y_up{3}(:,:), re_x_down{3}(:,:), re_y_down{3}(:,:), ...
% %                           re_x_left{3}(:,:), re_y_left{3}(:,:),re_x_right{3}(:,:),re_y_right{3}(:,:) );
% %                       
%                       
% % outfile = [BASE_PATH 'ProtocolSpecific\MOOG\RVOR_Pursuit\Que_Pursuit_World.dat'];
% % outfile = [BASE_PATH 'ProtocolSpecific\MOOG\3Dtuning\AzraelEye_Tra_Vet.dat'];
% % outfile = [BASE_PATH 'ProtocolSpecific\MOOG\3Dtuning\AzraelEye_Tra_Vis.dat'];
% 
% % outfile = [BASE_PATH 'ProtocolSpecific\MOOG\3Dtuning\QueEye_Tra_Vis.dat'];
% % outfile = [BASE_PATH 'ProtocolSpecific\MOOG\3Dtuning\ZebulonEye_Tra_Vis.dat'];
% 
% % outfile = [BASE_PATH 'ProtocolSpecific\MOOG\3Dtuning\QueLaby_Tra_Vis.dat'];
% outfile = [BASE_PATH 'ProtocolSpecific\MOOG\3Dtuning\ZebulonLaby_Tra_Vis.dat'];
% 
% printflag = 0;
% if (exist(outfile, 'file') == 0)    %file does not yet exist
%     printflag = 1;
% end
% fid = fopen(outfile, 'a');
% if (printflag)
%     fprintf(fid, 'FILE\t');
%     fprintf(fid, '\r\n');
% end
% fprintf(fid, '%s', buff);
% fprintf(fid, '\r\n');
% fclose(fid);
% 
% 
% 
% return;

function DirectionTuning3D_eyetrace(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

TEMPO_Defs;
Path_Defs;
ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for figure writing , if running batch, comment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % From here RVOR/pursuit
% if size(data.eye_data,1)>6
%     LEFT_EYE_1_2=9;
%     RIGHT_EYE_3_4=10;
% else
%     LEFT_EYE_1_2=7;
%     RIGHT_EYE_3_4=8;
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %************************************************************************%
% %plot Vertical vs. Horizontal
% switch (data.eye_flag)
%     case (LEFT_EYE_1_2)
%         Eye_Select='Left Eye'
%         Hor=1;        Ver=2;
%     case(RIGHT_EYE_3_4)
%         Eye_Select='Right Eye'
%         Hor=3;        Ver=4;
% end
% % 

temp_azimuth = data.moog_params(AZIMUTH,:,MOOG);
temp_elevation = data.moog_params(ELEVATION,:,MOOG);
temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG);
temp_spike_rates = data.spike_rates(SpikeChan, :); 
temp_total_trials = data.misc_params(OUTCOME, :);
temp_spike_data = data.spike_data(1,:);   % spike rasters
% temp_fp_rotate = data.moog_params(FP_ROTATE,:,MOOG);

% null_trials = logical( (temp_azimuth == data.one_time_params(NULL_VALUE)) );
null_trials = logical( (temp_elevation == data.one_time_params(NULL_VALUE)) );
%now, remove trials from direction and spike_rates that do not fall between BegTrial and EndTrial
% trials = 1:length(temp_azimuth);
trials = 1:length(temp_elevation);
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );
azimuth = temp_azimuth(~null_trials & select_trials);
elevation = temp_elevation(~null_trials & select_trials);
stim_type = temp_stim_type(~null_trials & select_trials);
spike_rates = temp_spike_rates(~null_trials & select_trials);
% fp_rotate = temp_fp_rotate(~null_trials & select_trials);

unique_azimuth = munique(azimuth');
unique_elevation = munique(elevation');
unique_stim_type = munique(stim_type');
% unique_fp_rotate = munique(fp_rotate');

% repeat = floor( length(temp_spike_rates) / (length(unique_stim_type)*(length(unique_fp_rotate)*(length(unique_elevation)-2)+2)+1) );
%Yong's bad Aihua's good for calculate repeat

trials_per_rep = (length(unique_azimuth)*length(unique_elevation)-14) * length(unique_stim_type) + 1;
repeat = floor( (EndTrial-(BegTrial-1)) / trials_per_rep);

for k = 1 : length(unique_stim_type)
    for i = 1 : length(unique_azimuth)
        for j = 1 : length(unique_elevation)
            select = find( temp_azimuth==unique_azimuth(i) & temp_elevation==unique_elevation(j) & temp_stim_type==unique_stim_type(k) );
            select_up = find( temp_azimuth==0 & temp_elevation==-90 & temp_stim_type==1 );
            if sum(select)>0
                for jj = 1 : length(select) % for convenience with sacrefice of some of the trials
%                Right Eye Que_Rotation, old_Que_c1_c71(Tra), Zeblon,
%                Lothar, Que_c290_c333(Tra), Que_Labyrinthectomy,
%                Zebulon_Laby
                    offset_x = mean( data.eye_data(3,201:300,select(jj)) ); % horizontal
                    offset_y = mean( data.eye_data(4,201:300,select(jj)) ); % vertical
                    resp_x{k,i,j}(jj,:) = data.eye_data(3,201:600,select(jj)) - offset_x;  % horizontal   
                    resp_y{k,i,j}(jj,:) = data.eye_data(4,201:600,select(jj)) - offset_y;  % horizontal 
                    resp_yaw1L{k,i,j}(jj,:) = data.eye_data(3,324:523,select(jj)) - 0;  % horizontal   
                    resp_yaw2L{k,i,j}(jj,:) = data.eye_data(4,324:523,select(jj)) - 0;  % horizontal.
%                Left Eye Azrael, Lothar(after..),
%                new_Que_c86__butCh.3,4 are O.K., too
%                     offset_x = mean( data.eye_data(1,201:300,select(jj)) ); % horizontal
%                     offset_y = mean( data.eye_data(2,201:300,select(jj)) ); % vertical
%                     resp_x{k,i,j}(jj,:) = data.eye_data(1,201:600,select(jj)) - offset_x;  % horizontal   
%                     resp_y{k,i,j}(jj,:) = data.eye_data(2,201:600,select(jj)) - offset_y;  % horizontal
%                Eye select (for figure wrinting)
%                     offset_x = mean( data.eye_data(Hor,201:300,select(jj)) ); % horizontal
%                     offset_y = mean( data.eye_data(Ver,201:300,select(jj)) ); % vertical
%                     resp_x{k,i,j}(jj,:) = data.eye_data(Hor,201:600,select(jj)) - offset_x;  % horizontal   
%                     resp_y{k,i,j}(jj,:) = data.eye_data(Ver,201:600,select(jj)) - offset_y;  % vertical
                end
            else
                resp_x{k,i,j}(:,:) = resp_x{k,1,j}(:,:);%??? +-90 elevation
                resp_y{k,i,j}(:,:) = resp_y{k,1,j}(:,:);
            end
        end
    end
    
    % resp_x{k,i,1}(:,:) i=1 means stim_type = vestibular, if you want to
    % visual turn to resp_x{k,2,1}(:,:)
    
    resp_x_up{k}(:,:) = resp_x{k,1,1}(:,:);     resp_y_up{k}(:,:) = resp_y{k,1,1}(:,:);
    resp_x_down{k}(:,:) = resp_x{k,1,5}(:,:);   resp_y_down{k}(:,:) = resp_y{k,1,5}(:,:);
    resp_x_left{k}(:,:) = resp_x{k,5,3}(:,:);   resp_y_left{k}(:,:) = resp_y{k,5,3}(:,:);
    resp_x_right{k}(:,:) = resp_x{k,1,3}(:,:);  resp_y_right{k}(:,:) = resp_y{k,1,3}(:,:);  
%     resp_x_up{k}(:,:) = resp_x{k,2,5}(:,:);     resp_y_up{k}(:,:) = resp_y{k,2,5}(:,:);
%     resp_x_down{k}(:,:) = resp_x{k,2,1}(:,:);   resp_y_down{k}(:,:) = resp_y{k,2,1}(:,:);
%     resp_x_left{k}(:,:) = resp_x{k,2,3}(:,:);   resp_y_left{k}(:,:) = resp_y{k,2,3}(:,:);
%     resp_x_right{k}(:,:) = resp_x{k,2,7}(:,:);  resp_y_right{k}(:,:) = resp_y{k,2,7}(:,:);  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%   Yong's Eye analysis matlab file  No more use
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % reshape the matrix for output text file
    %% 400points-400-400-400...400*5repeats=2000points only for fixation
    
    
%     re_x_up{k} = reshape(resp_x_up{k}(:,:)',1,repeat*400);
%     re_y_up{k} = reshape(resp_y_up{k}(:,:)',1,repeat*400);
%     re_x_down{k} = reshape(resp_x_down{k}(:,:)',1,repeat*400);
%     re_y_down{k} = reshape(resp_y_down{k}(:,:)',1,repeat*400);
%     re_x_left{k} = reshape(resp_x_left{k}(:,:)',1,repeat*400);
%     re_y_left{k} = reshape(resp_y_left{k}(:,:)',1,repeat*400);
%     re_x_right{k} = reshape(resp_x_right{k}(:,:)',1,repeat*400);
%     re_y_right{k} = reshape(resp_y_right{k}(:,:)',1,repeat*400);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     mean for 5 repeatition for figures
    re_x_up{k} = mean(resp_x_up{k}(:,:));
    re_y_up{k} = mean(resp_y_up{k}(:,:));
    re_x_down{k} = mean(resp_x_down{k}(:,:));
    re_y_down{k} = mean(resp_y_down{k}(:,:));
    re_x_left{k} = mean(resp_x_left{k}(:,:));
    re_y_left{k} = mean(resp_y_left{k}(:,:));
    re_x_right{k} = mean(resp_x_right{k}(:,:));
    re_y_right{k} = mean(resp_y_right{k}(:,:));
%     

end
% % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% %For figure writing, in roder to batch comment!!!!
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% 
% 
%    title1 = 'up';
%     title2 = 'down';
%     title3 = 'left';
%     title4 = 'right';
%     
%     
%     figure(2)
% 
% subplot(4,1,1)
% plot(re_x_up{1},'r.');
%     hold on;
%     xlim( [1, 400] );
%     set(gca, 'XTickMode','manual');
%     set(gca, 'xtick',[1,100,200,300,400]);
%     set(gca, 'xticklabel','0|0.5|1|1.5|2'); 
% %     ylim([-20, 20]);
%     ylabel('(deg)');
%     title(['Eye Position /  ',title1]);
% 
%     plot(re_y_up{1},'b.');
%     hold off;
% 
% 
% subplot(4,1,2)
% plot(re_x_down{1},'r.');
%     hold on;
%     xlim( [1, 400] );
%     set(gca, 'XTickMode','manual');
%     set(gca, 'xtick',[1,100,200,300,400]);
%     set(gca, 'xticklabel','0|0.5|1|1.5|2'); 
% %     ylim([-20, 20]);
%     ylabel('(deg)');
%     title(['Eye Position /  ',title2]);
% 
%     plot(re_y_down{1},'b.');
%     hold off;
% 
% subplot(4,1,3)
% plot(re_x_left{1},'r.');
%     hold on;
%     xlim( [1, 400] );
%     set(gca, 'XTickMode','manual');
%     set(gca, 'xtick',[1,100,200,300,400]);
%     set(gca, 'xticklabel','0|0.5|1|1.5|2'); 
% %     ylim([-20, 20]);
%     ylabel('(deg)');
%     title(['Eye Position /  ',title3]);
% 
%     plot(re_y_left{1},'b.');
%     hold off;
%     
% subplot(4,1,4)
% plot(re_x_right{1},'r.');
%     hold on;
%     xlim( [1, 400] );
%     set(gca, 'XTickMode','manual');
%     set(gca, 'xtick',[1,100,200,300,400]);
%     set(gca, 'xticklabel','0|0.5|1|1.5|2'); 
% %     ylim([-20, 20]);
%     ylabel('(deg)');
%     title(['Eye Position /  ',title4]);
% 
%     plot(re_y_right{1},'b.');
%     hold off;
%     
%  figure(3)
% 
% subplot(4,1,1)
% plot(re_x_up{2},'r.');
%     hold on;
%     xlim( [1, 400] );
%     set(gca, 'XTickMode','manual');
%     set(gca, 'xtick',[1,100,200,300,400]);
%     set(gca, 'xticklabel','0|0.5|1|1.5|2'); 
% %     ylim([-20, 20]);
%     ylabel('(deg)');
%     title(['Eye Position /  ',title1]);
% 
%     plot(re_y_up{2},'b.');
%     hold off;
% 
% 
% subplot(4,1,2)
% plot(re_x_down{2},'r.');
%     hold on;
%     xlim( [1, 400] );
%     set(gca, 'XTickMode','manual');
%     set(gca, 'xtick',[1,100,200,300,400]);
%     set(gca, 'xticklabel','0|0.5|1|1.5|2'); 
% %     ylim([-20, 20]);
%     ylabel('(deg)');
%     title(['Eye Position /  ',title2]);
% 
%     plot(re_y_down{2},'b.');
%     hold off;
% 
% subplot(4,1,3)
% plot(re_x_left{2},'r.');
%     hold on;
%     xlim( [1, 400] );
%     set(gca, 'XTickMode','manual');
%     set(gca, 'xtick',[1,100,200,300,400]);
%     set(gca, 'xticklabel','0|0.5|1|1.5|2'); 
% %     ylim([-20, 20]);
%     ylabel('(deg)');
%     title(['Eye Position /  ',title3]);
% 
%     plot(re_y_left{2},'b.');
%     hold off;
%     
% subplot(4,1,4)
% plot(re_x_right{2},'r.');
%     hold on;
%     xlim( [1, 400] );
%     set(gca, 'XTickMode','manual');
%     set(gca, 'xtick',[1,100,200,300,400]);
%     set(gca, 'xticklabel','0|0.5|1|1.5|2'); 
% %     ylim([-20, 20]);
%     ylabel('(deg)');
%     title(['Eye Position /  ',title4]);
% 
%     plot(re_y_right{2},'b.');
%     hold off;
%     
%     
%     
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                         output to text file
sprint_txt = ['%s'];
for i = 1 : 400*repeat*10
    sprint_txt = [sprint_txt, ' %1.2f'];    
end


% if you want to save 'world', select 2, 'pursuit', select 3


% buff= sprintf(sprint_txt, FILE, repeat, re_y_up{1}(:,:), re_y_down{1}(:,:), ...
%                           re_x_left{1}(:,:), re_x_right{1}(:,:) );
% % buff= sprintf(sprint_txt, FILE, repeat, re_x_up{2}(:,:), re_y_up{2}(:,:), re_x_down{2}(:,:), re_y_down{2}(:,:), ...
% %                           re_x_left{2}(:,:), re_y_left{2}(:,:),re_x_right{2}(:,:),re_y_right{2}(:,:) );
% % buff= sprintf(sprint_txt, FILE, repeat, re_x_up{3}(:,:), re_y_up{3}(:,:), re_x_down{3}(:,:), re_y_down{3}(:,:), ...
% %                           re_x_left{3}(:,:), re_y_left{3}(:,:),re_x_right{3}(:,:),re_y_right{3}(:,:) );
% %                       
%                       
% % outfile = [BASE_PATH 'ProtocolSpecific\MOOG\RVOR_Pursuit\Que_Pursuit_World.dat'];
% % outfile = [BASE_PATH 'ProtocolSpecific\MOOG\3Dtuning\AzraelEye_Tra_Vet.dat'];
% % outfile = [BASE_PATH 'ProtocolSpecific\MOOG\3Dtuning\AzraelEye_Tra_Vis.dat'];
% 
% % outfile = [BASE_PATH 'ProtocolSpecific\MOOG\3Dtuning\QueEye_Tra_Vet.dat'];
% % outfile = [BASE_PATH 'ProtocolSpecific\MOOG\3Dtuning\ZebulonEye_Tra_Vet.dat'];
% 
% % outfile = [BASE_PATH 'ProtocolSpecific\MOOG\3Dtuning\QueLaby_Tra_Vet.dat'];
% outfile = [BASE_PATH 'ProtocolSpecific\MOOG\3Dtuning\ZebulonLaby_Tra_Vet.dat'];
buff= sprintf(sprint_txt, FILE, resp_yaw1L{1,1,1}(:,:), resp_yaw2L{1,1,1}(:,:), resp_yaw1L{1,1,5}(:,:), resp_yaw2L{1,1,5}(:,:));
outfile = ['Z:\Users\Yong\MSTrotationdarkness_eyeposition.dat'];

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



% %%%%%                         output to text file
% sprint_txt2 = ['%s'];
% for i = 1 : 400*repeat*10
%     sprint_txt2 = [sprint_txt2, ' %1.2f'];    
% end
% 
% 
% % if you want to save 'world', select 2, 'pursuit', select 3
% 
% 
% % buff= sprintf(sprint_txt2, FILE, repeat, re_x_up{1}(:,:), re_y_up{1}(:,:), re_x_down{1}(:,:), re_y_down{1}(:,:), ...
% %                           re_x_left{1}(:,:), re_y_left{1}(:,:),re_x_right{1}(:,:),re_y_right{1}(:,:) );
% buff= sprintf(sprint_txt, FILE, repeat, re_y_up{2}(:,:), re_y_down{2}(:,:), ...
%                           re_x_left{2}(:,:), re_x_right{2}(:,:) );
% % buff= sprintf(sprint_txt, FILE, repeat, re_x_up{3}(:,:), re_y_up{3}(:,:), re_x_down{3}(:,:), re_y_down{3}(:,:), ...
% %                           re_x_left{3}(:,:), re_y_left{3}(:,:),re_x_right{3}(:,:),re_y_right{3}(:,:) );
% %                       
%                       
% % outfile = [BASE_PATH 'ProtocolSpecific\MOOG\RVOR_Pursuit\Que_Pursuit_World.dat'];
% % outfile = [BASE_PATH 'ProtocolSpecific\MOOG\3Dtuning\AzraelEye_Tra_Vet.dat'];
% % outfile = [BASE_PATH 'ProtocolSpecific\MOOG\3Dtuning\AzraelEye_Tra_Vis.dat'];
% 
% % outfile = [BASE_PATH 'ProtocolSpecific\MOOG\3Dtuning\QueEye_Tra_Vis.dat'];
% % outfile = [BASE_PATH 'ProtocolSpecific\MOOG\3Dtuning\ZebulonEye_Tra_Vis.dat'];
% 
% % outfile = [BASE_PATH 'ProtocolSpecific\MOOG\3Dtuning\QueLaby_Tra_Vis.dat'];
% outfile = [BASE_PATH 'ProtocolSpecific\MOOG\3Dtuning\ZebulonLaby_Tra_Vis.dat'];
% 
% printflag = 0;
% if (exist(outfile, 'file') == 0)    %file does not yet exist
%     printflag = 1;
% end
% fid = fopen(outfile, 'a');
% if (printflag)
%     fprintf(fid, 'FILE\t');
%     fprintf(fid, '\r\n');
% end
% fprintf(fid, '%s', buff);
% fprintf(fid, '\r\n'); 
% fclose(fid);



return;



