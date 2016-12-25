function Rotation3D_eyetrace(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE); 

TEMPO_Defs;
Path_Defs;
ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% From here RVOR/pursuit
if size(data.eye_data,1)>6
    LEFT_EYE_1_2=9;
    RIGHT_EYE_3_4=10;
else
    LEFT_EYE_1_2=7;
    RIGHT_EYE_3_4=8;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%************************************************************************%
%plot Vertical vs. Horizontal
switch (data.eye_flag)
    case (LEFT_EYE_1_2)
        Eye_Select='Left Eye'
        Hor=1;        Ver=2;
    case(RIGHT_EYE_3_4)
        Eye_Select='Right Eye'
        Hor=3;        Ver=4;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%get the Sample Rate
filename = [PATH FILE];
fid = htbOpen(filename);
h = htbGetHd(fid, 1); %Make Sure: Database#1: Eye Traces
SR=(h.speed_units/h.speed)/(h.skip+1);%Sample Rate
clear filename fid h ;
%%%%%%%%%%%%%%%%%%%%%%%%%%

StartPoint=SR*1+1;EndPoint=SR*3;
NewX(:,:)=data.eye_data(Hor,StartPoint:EndPoint,:);
NewY(:,:)=data.eye_data(Ver,StartPoint:EndPoint,:);
% figure(2);plot(NewX, NewY,'b');xlabel('Horizontal Position');ylabel('Vertical Position');title([FILE,'/',Eye_Select]);

%************************************************************************%
%position vs. time and velocity vs. time according to different direction
TimeLabel=(StartPoint:EndPoint)/SR*1000;%ms
temp_rot_elevation=data.moog_params(ROT_ELEVATION,:,MOOG);
temp_fp_rotate = data.moog_params(FP_ROTATE,:,MOOG);
%Add by Katsu 11/03/06
temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG);
% for k=1:length(temp_stim_type)
for i=1:8
   clear Select_Trial; Select_Trial=find(temp_fp_rotate==1 & temp_rot_elevation==(i-1)*45);%FP_ROTATE=1 means world-fixed???
   %FP_rotate=0 means head-fixed now select 0
    clear Position_H1 Position_H2  Position_TH1 Position_TH2 Position_V1 Position_V2  Position_TV1 Position_TV2;
    Position_H1(:,:)=data.eye_data(Hor,StartPoint:EndPoint,Select_Trial);%Horizontal_Rotation
    Position_TH1(:,:)=data.eye_data(5,StartPoint:EndPoint,Select_Trial);%Target_Horizontal_Rotation
    Position_V1(:,:)=data.eye_data(Ver,StartPoint:EndPoint,Select_Trial);%Vertical_Rotation   
    Position_TV1(:,:)=data.eye_data(6,StartPoint:EndPoint,Select_Trial);%Target_Vertical_Rotation
  
    clear Select_Trial;Select_Trial=find(temp_fp_rotate==2 & temp_rot_elevation==(i-1)*45);%FP_ROTATE=2 means pursuit only???
    Position_H2(:,:)=data.eye_data(Hor,StartPoint:EndPoint, Select_Trial);%Horizontal_Pursuit
    Position_TH2(:,:)=data.eye_data(5,StartPoint:EndPoint,Select_Trial);%Target_Horizontal_Pursuit
    Position_V2(:,:)=data.eye_data(Ver,StartPoint:EndPoint,Select_Trial);%Vertical_Pursuit
    Position_TV2(:,:)=data.eye_data(6,StartPoint:EndPoint,Select_Trial);%Target_Vertical_Pursuit   

    figure(3);
    set(3, 'Name', 'World-fix');
    subplot(2,4,i);plot(TimeLabel,Position_H1(:,1), 'b', TimeLabel,Position_H2(:,1),'r',TimeLabel,Position_TH1(:,1),'k',TimeLabel,Position_TH2,'c'); 
    %legend('Rotation','Pursuit only','Target // Rotation','Target // Pursuit only',2); 
%     hold on;plot(TimeLabel,Position_H1(:,2:end), 'b', TimeLabel,Position_H2(:,2:end),'r',TimeLabel,Position_TH1(:,2:end),'k',TimeLabel,Position_TH2(:,2:end),'c'); 
    xlabel('Time (ms)');ylabel('Horizontal Position');title([FILE,'/',Eye_Select,'/',num2str((i-1)*45)]);axis([1000 3000 -12 12]);
 
    figure(2);
    set(2, 'Name', 'World-fix');
    subplot(2,4,i);plot(TimeLabel,Position_H1(:,:), 'b', TimeLabel,Position_H2(:,:),'r',TimeLabel,Position_TH1(:,:),'k',TimeLabel,Position_TH2,'c'); 
    %legend('Rotation','Pursuit only','Target // Rotation','Target // Pursuit only',2); 
%     hold on;plot(TimeLabel,Position_H1(:,2:end), 'b', TimeLabel,Position_H2(:,2:end),'r',TimeLabel,Position_TH1(:,2:end),'k',TimeLabel,Position_TH2(:,2:end),'c'); 
    xlabel('Time (ms)');ylabel('Horizontal Position');title([FILE,'/',Eye_Select,'/',num2str((i-1)*45)]);axis([1000 3000 -12 12]);
    
figure(4);subplot(2,4,i);plot(TimeLabel,Position_V1(:,:), 'b', TimeLabel,Position_V2(:,:),'r',TimeLabel,Position_TV1(:,:),'k',TimeLabel,Position_TV2(:,:),'c'); 
    %legend('Rotation','Pursuit only','Target // Rotation','Target // Pursuit only',2); 
%     hold on;plot(TimeLabel,Position_V1(:,2:end), 'b', TimeLabel,Position_V2(:,2:end),'r',TimeLabel,Position_TV1(:,2:end),'k',TimeLabel,Position_TV2(:,2:end),'c'); 
    xlabel('Time (ms)');ylabel('Vertical Position');title([FILE,'/',Eye_Select,'/',num2str((i-1)*45)]);axis([1000 3000 -12 12]);    

    
    clear Velocity_H1 Velocity_H2  Velocity_TH1 Velocity_TH2 Velocity_V1 Velocity_V2 Velocity_TV1 Velocity_TV2 ;      
    for j=1:size(Position_H1,2)
        Velocity_H1(:,j)=fderiv(Position_H1(:,j),15,SR);%Horizontal_Rotation
        Velocity_H2(:,j)=fderiv(Position_H2(:,j),15,SR);%Horizontal_Pursuit
        Velocity_TH1(:,j)=fderiv(Position_TH1(:,j),15,SR);%Target_Horizontal_Rotation
        Velocity_TH2(:,j)=fderiv(Position_TH2(:,j),15,SR);%Target_Horizontal_Pursuit
        
        Velocity_V1(:,j)=fderiv(Position_V1(:,j),15,SR);%Vertical_Rotation
        Velocity_V2(:,j)=fderiv(Position_V2(:,j),15,SR);%Vertical_Pursuit
        Velocity_TV1(:,j)=fderiv(Position_TV1(:,j),15,SR);%Target_Vertical_Rotation
        Velocity_TV2(:,j)=fderiv(Position_TV2(:,j),15,SR);%Target_Vertical_Pursuit        
    end    
    figure(5);subplot(2,4,i);plot(TimeLabel,Velocity_H1(:,:), 'b', TimeLabel,Velocity_H2(:,:),'r',TimeLabel,Velocity_TH1(:,:),'k',TimeLabel, Velocity_TH2(:,:),'c');
    %legend('Rotation','Pursuit only','Target // Rotation','Target // Pursuit only',2); 
    hold on; %plot(TimeLabel,Velocity_H1(:,2:end), 'b', TimeLabel,Velocity_H2(:,2:end),'r',TimeLabel,Velocity_TH1(:,2:end),'k',TimeLabel, Velocity_TH2(:,2:end),'c');
    xlabel('Time (ms)');ylabel(' Horizontal Velocity');title([FILE,'/',Eye_Select,'/',num2str((i-1)*45)]);axis([1000 3000 -30 30]);   
    
    figure(6);subplot(2,4,i);plot(TimeLabel,Velocity_V1(:,1), 'b', TimeLabel,Velocity_V2(:,1),'r',TimeLabel,Velocity_TV1(:,1),'k',TimeLabel, Velocity_TV2,'c');
    %legend('Rotation','Pursuit only','Target // Rotation','Target // Pursuit only',2); 
    hold on; plot(TimeLabel,Velocity_V1(:,2:end), 'b', TimeLabel,Velocity_V2(:,2:end),'r',TimeLabel,Velocity_TV1(:,2:end),'k',TimeLabel, Velocity_TV2(:,2:end),'c');
    xlabel('Time (ms)');ylabel(' Vertical Velocity');title([FILE,'/',Eye_Select,'/',num2str((i-1)*45)]);axis([1000 3000 -30 30]);       
end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   above are made by Aihua for RVOR/pursuit
%   Following is by Yong for Eye Trace in Rotation and Translation
% 
% % the following for Katzu, temporarilly
% temp_azimuth = data.moog_params(ROT_AZIMUTH,:,MOOG);
% temp_elevation = data.moog_params(ROT_ELEVATION,:,MOOG);
% temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG);
% temp_spike_rates = data.spike_rates(SpikeChan, :); 
% temp_total_trials = data.misc_params(OUTCOME, :);
% temp_spike_data = data.spike_data(1,:);   % spike rasters
% 
% null_trials = logical( (temp_azimuth == data.one_time_params(NULL_VALUE)) );
% %now, remove trials from direction and spike_rates that do not fall between BegTrial and EndTrial
% trials = 1:length(temp_azimuth);
% select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );
% azimuth = temp_azimuth(~null_trials & select_trials);
% elevation = temp_elevation(~null_trials & select_trials);
% stim_type = temp_stim_type(~null_trials & select_trials);
% spike_rates = temp_spike_rates(~null_trials & select_trials);
% 
% unique_azimuth = munique(azimuth');
% unique_elevation = munique(elevation');
% unique_stim_type = munique(stim_type');
% 
% repeat = floor( length(temp_spike_rates) / (length(unique_stim_type)*(length(unique_azimuth)*(length(unique_elevation)-2)+2)+1) );
% for k = 1 : length(unique_stim_type)
%     for i = 1 : length(unique_azimuth)
%         for j = 1 : length(unique_elevation)
%             select = find( temp_azimuth==unique_azimuth(i) & temp_elevation==unique_elevation(j) & temp_stim_type==unique_stim_type(k) );
%             if sum(select)>0
%                 for jj = 1 : repeat % for convenience with sacrefice of some of the trials
%                     offset_x = mean( data.eye_data(3,201:300,select(jj)) ); % horizontal
%                     offset_y = mean( data.eye_data(4,201:300,select(jj)) ); % vertical
%                     resp_x{k,i,j}(jj,:) = data.eye_data(3,201:600,select(jj)) - offset_x;  % horizontal   
%                     resp_y{k,i,j}(jj,:) = data.eye_data(4,201:600,select(jj)) - offset_y;  % horizontal   
%                 end
%             else
%                 resp_x{k,i,j}(:,:) = resp_x{k,1,j}(:,:);
%                 resp_y{k,i,j}(:,:) = resp_y{k,1,j}(:,:);
%             end
%         end
%     end
%     resp_x_up{k}(:,:) = resp_x{k,1,1}(:,:);     resp_y_up{k}(:,:) = resp_y{k,1,1}(:,:);
%     resp_x_down{k}(:,:) = resp_x{k,1,5}(:,:);   resp_y_down{k}(:,:) = resp_y{k,1,5}(:,:);
%     resp_x_left{k}(:,:) = resp_x{k,5,3}(:,:);   resp_y_left{k}(:,:) = resp_y{k,5,3}(:,:);
%     resp_x_right{k}(:,:) = resp_x{k,1,3}(:,:);  resp_y_right{k}(:,:) = resp_y{k,1,3}(:,:);  
%     % reshape the matrix for output text file
%     re_x_up{k} = reshape(resp_x_up{k}(:,:)',1,repeat*400);
%     re_y_up{k} = reshape(resp_y_up{k}(:,:)',1,repeat*400);
%     re_x_down{k} = reshape(resp_x_down{k}(:,:)',1,repeat*400);
%     re_y_down{k} = reshape(resp_y_down{k}(:,:)',1,repeat*400);
%     re_x_left{k} = reshape(resp_x_left{k}(:,:)',1,repeat*400);
%     re_y_left{k} = reshape(resp_y_left{k}(:,:)',1,repeat*400);
%     re_x_right{k} = reshape(resp_x_right{k}(:,:)',1,repeat*400);
%     re_y_right{k} = reshape(resp_y_right{k}(:,:)',1,repeat*400);
% end
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % output to text file
% sprint_txt = ['%s'];
% for i = 1 : 400*repeat*10
%     sprint_txt = [sprint_txt, ' %1.2f'];    
% end
% % if you want to save 'vestibular', select 1, 'visual', select 2
% % buff= sprintf(sprint_txt, FILE, repeat, re_x_up{1}(:,:), re_y_up{1}(:,:), re_x_down{1}(:,:), re_y_down{1}(:,:), ...
% %                           re_x_left{1}(:,:), re_y_left{1}(:,:),re_x_right{1}(:,:),re_y_right{1}(:,:) );
% buff= sprintf(sprint_txt, FILE, repeat, re_x_up{2}(:,:), re_y_up{2}(:,:), re_x_down{2}(:,:), re_y_down{2}(:,:), ...
%                           re_x_left{2}(:,:), re_y_left{2}(:,:),re_x_right{2}(:,:),re_y_right{2}(:,:) );
% % outfile = [BASE_PATH 'ProtocolSpecific\MOOG\Rotation3D\Eye_rot_ves.dat'];
% outfile = [BASE_PATH 'ProtocolSpecific\MOOG\Rotation3D\Eye_rot_vis.dat'];
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

