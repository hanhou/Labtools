
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%       11/28/06        accelerometer  CH. 5 and 6
%%     Hor=5; Ver=6;
%%%            Right eye should be selected

function Rotation3D_eyetrace(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE); 

TEMPO_Defs;
Path_Defs;
ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP

% FOr Que 
% LEFT_EYE_1_2=7;
% RIGHT_EYE_3_4=8;
%For Lothar
LEFT_EYE_1_2=9;
RIGHT_EYE_3_4=10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%************************************************************************%
%plot Vertical vs. Horizontal
switch (data.eye_flag)
    case (LEFT_EYE_1_2)
        Eye_Select='Left Eye'
        Hor=1;        Ver=2;
    case(RIGHT_EYE_3_4)
        Eye_Select='Right Eye'
        Hor=5;        Ver=6;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%get the Sample Rate
filename = [PATH FILE];
fid = htbOpen(filename);
h = htbGetHd(fid, 1); %Make Sure: Database#1: Eye Traces
SR=(h.speed_units/h.speed)/(h.skip+1)%Sample Rate
clear filename fid h ;
%%%%%%%%%%%%%%%%%%%%%%%%%%

% StartPoint=SR*1+1;EndPoint=SR*3;
StartPoint=1;EndPoint=999;

%position vs. time and velocity vs. time according to different direction
TimeLabel=(StartPoint:EndPoint)/SR*1000;%ms
% temp_rot_elevation=data.moog_params(ROT_ELEVATION,:,MOOG);
% temp_fp_rotate = data.moog_params(FP_ROTATE,:,MOOG);
temp_rot_elevation=data.moog_params(ROT_AZIMUTH,:,MOOG);%That's tric! only change ELE to AZI
temp_fp_rotate = data.moog_params(TT_MODE,:,MOOG);%Same Tric! only change FP-ROTATE to TT_MODE
temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG);

% clear Select_Trial; Select_Trial=find(temp_stim_type==1 & temp_fp_rotate==0 & temp_rot_elevation==0); %0 degrees representative
% for Que
%     Position_TH1(:,:)=data.eye_data(5,StartPoint:EndPoint,Select_Trial);%Target_Horizontal_Rotation
%     Position_TV1(:,:)=data.eye_data(6,StartPoint:EndPoint,Select_Trial);%Target_Vertical_Rotation
% for Lothar    
%     Position_TH1(:,:)=data.eye_data(7,StartPoint:EndPoint,Select_Trial);%Target_Horizontal_Rotation
%     Position_TV1(:,:)=data.eye_data(8,StartPoint:EndPoint,Select_Trial);%Target_Vertical_Rotation
    %%%%%%% Position data is units, which are from tempo side, it should be
    %%%%%%% converted.
%         Position_TH1(:,:)=Position_TH1(:,:)/655.36;%+or-5V; 2^15=32768units; corresponds to whole display 100degrees/2 (because+-5V)
%         Position_TV1(:,:)=Position_TV1(:,:)/655.36;%so 1 degree=655.36 units
        
%%% This version only pick-up 0 deg and 180 deg pitch-up and pitch-down

for i=1:2
    clear Position_H1 Position_H2  Position_H3 Position_H4 Position_V1 Position_V2 Position_V3 Position_V4  %Position_H0  Position_V0;
   
    clear Select_Trial; Select_Trial=find(temp_stim_type==1 & temp_fp_rotate==0 & temp_rot_elevation==(i-1)*180);%FP_ROTATE(TT_MODA)=0 Tilt=Trans 
    Position_H1(:,:)=data.eye_data(Hor,StartPoint:EndPoint,Select_Trial);%Horizontal
    Position_V1(:,:)=data.eye_data(Ver,StartPoint:EndPoint,Select_Trial);%Vertical
    Position_H1=mean(Position_H1');
    Position_V1=mean(Position_V1');
  
    clear Select_Trial; Select_Trial=find(temp_stim_type==1 & temp_fp_rotate==1 & temp_rot_elevation==(i-1)*180);%FP_ROTATE=1 Tilt-Trans
    Position_H2(:,:)=data.eye_data(Hor,StartPoint:EndPoint, Select_Trial);%Horizontal_
    Position_V2(:,:)=data.eye_data(Ver,StartPoint:EndPoint,Select_Trial);%Vertical_
    Position_H2=mean(Position_H2');
    Position_V2=mean(Position_V2');
    
    clear Select_Trial; Select_Trial=find(temp_stim_type==1 & temp_fp_rotate==2 & temp_rot_elevation==(i-1)*180);%FP_ROTATE=2 Tilt only
    Position_H3(:,:)=data.eye_data(Hor,StartPoint:EndPoint, Select_Trial);%Horizontal_
    Position_V3(:,:)=data.eye_data(Ver,StartPoint:EndPoint,Select_Trial);%Vertical_
    Position_H3=mean(Position_H3');
    Position_V3=mean(Position_V3');
    
    clear Select_Trial; Select_Trial=find(temp_stim_type==1 & temp_fp_rotate==3 & temp_rot_elevation==(i-1)*180);%FP_ROTATE=3 Trans only
    Position_H4(:,:)=data.eye_data(Hor,StartPoint:EndPoint, Select_Trial);%Horizontal_
    Position_V4(:,:)=data.eye_data(Ver,StartPoint:EndPoint,Select_Trial);%Vertical_
    Position_H4=mean(Position_H4');
    Position_V4=mean(Position_V4');
 


    figure(2);set(2,'Position', [5,5 980,650], 'name',FILE);orient landscape;% middle 2 sec
    subplot(4,2,i);plot(TimeLabel,Position_H1(:,:), 'r');
    title(['X/  Ro ',num2str((i-1)*180),'  (Pitch)  Tra  ',num2str((i-1)*180+90),' deg / Tilt+Trans']);set(gca, 'xtick', [] );axis([1000 3000 13 18]);
    hold on;
    subplot(4,2,i+2);plot(TimeLabel,Position_H2(:,:), 'r');set(gca, 'xtick', [] );title(['Tilt-Trans']);axis([1000 3000 13 18]);hold on;
    subplot(4,2,i+4);plot(TimeLabel,Position_H3(:,:), 'r');set(gca, 'xtick', [] );title(['Tilt']);axis([1000 3000 13 18]);hold on;
    subplot(4,2,i+6);plot(TimeLabel,Position_H4(:,:), 'r');set(gca, 'xtick', [] );title(['Trans']);axis([1000 3000 13 18]);hold on;
    
    figure(3);set(3,'Position', [5,5 980,650], 'name',FILE);orient landscape;
    subplot(4,2,i);plot(TimeLabel,Position_V1(:,:), 'b');
    title(['Y/  Ro ',num2str((i-1)*180),'  (Pitch)   Tra  ',num2str((i-1)*180+90),' deg / Tilt+Trans']);set(gca, 'xtick', [] );axis([1000 3000 13 18]);
    hold on;
    subplot(4,2,i+2);plot(TimeLabel,Position_V2(:,:), 'b');set(gca, 'xtick', [] );title(['Tilt-Trans']);axis([1000 3000 13 18]);hold on;
    subplot(4,2,i+4);plot(TimeLabel,Position_V3(:,:), 'b');set(gca, 'xtick', [] );title(['Tilt']);axis([1000 3000 13 18]);hold on;
    subplot(4,2,i+6);plot(TimeLabel,Position_V4(:,:), 'b');set(gca, 'xtick', [] );title(['Trans']);axis([1000 3000 13 18]);hold on;
    
    
    figure(4);set(4,'Position', [5,5 980,650], 'name',FILE);orient landscape;% wholoe 5 sec
    subplot(4,2,i);plot(TimeLabel,Position_H1(:,:), 'r');
    title(['X/  Ro ',num2str((i-1)*180),'  (Pitch)  Tra  ',num2str((i-1)*180+90),' deg / Tilt+Trans']);set(gca, 'xtick', [] );
    hold on;
    subplot(4,2,i+2);plot(TimeLabel,Position_H2(:,:), 'r');set(gca, 'xtick', [] );title(['Tilt-Trans']);hold on;
    subplot(4,2,i+4);plot(TimeLabel,Position_H3(:,:), 'r');set(gca, 'xtick', [] );title(['Tilt']);hold on;
    subplot(4,2,i+6);plot(TimeLabel,Position_H4(:,:), 'r');set(gca, 'xtick', [] );title(['Trans']);hold on;
    
    figure(5);set(5,'Position', [5,5 980,650], 'name',FILE);orient landscape;
    subplot(4,2,i);plot(TimeLabel,Position_V1(:,:), 'b');
    title(['Y/  Ro ',num2str((i-1)*180),'  (Pitch)   Tra  ',num2str((i-1)*180+90),' deg / Tilt+Trans']);set(gca, 'xtick', [] );
    hold on;
    subplot(4,2,i+2);plot(TimeLabel,Position_V2(:,:), 'b');set(gca, 'xtick', [] );title(['Tilt-Trans']);hold on;
    subplot(4,2,i+4);plot(TimeLabel,Position_V3(:,:), 'b');set(gca, 'xtick', [] );title(['Tilt']);hold on;
    subplot(4,2,i+6);plot(TimeLabel,Position_V4(:,:), 'b');set(gca, 'xtick', [] );title(['Trans']);hold on;
    
    
end

%%% This version only pick-up 90 deg and 270 deg Roll-right and Roll-left
for i=1:2
    clear Position_H5 Position_H6  Position_H7 Position_H8 Position_V5 Position_V6 Position_V7 Position_V8  %Position_H0  Position_V0;
   
    clear Select_Trial; Select_Trial=find(temp_stim_type==1 & temp_fp_rotate==0 & temp_rot_elevation==(i-1)*180+90);%FP_ROTATE(TT_MODA)=0 Tilt=Trans 
    Position_H5(:,:)=data.eye_data(Hor,StartPoint:EndPoint,Select_Trial);%Horizontal
    Position_V5(:,:)=data.eye_data(Ver,StartPoint:EndPoint,Select_Trial);%Vertical
    Position_H5=mean(Position_H5');
    Position_V5=mean(Position_V5');
  
    clear Select_Trial; Select_Trial=find(temp_stim_type==1 & temp_fp_rotate==1 & temp_rot_elevation==(i-1)*180+90);%FP_ROTATE=1 Tilt-Trans
    Position_H6(:,:)=data.eye_data(Hor,StartPoint:EndPoint, Select_Trial);%Horizontal
    Position_V6(:,:)=data.eye_data(Ver,StartPoint:EndPoint,Select_Trial);%Vertical
    Position_H6=mean(Position_H6');
    Position_V6=mean(Position_V6');
    
    clear Select_Trial; Select_Trial=find(temp_stim_type==1 & temp_fp_rotate==2 & temp_rot_elevation==(i-1)*180+90);%FP_ROTATE=2 Tilt only
    Position_H7(:,:)=data.eye_data(Hor,StartPoint:EndPoint, Select_Trial);%Horizontal
    Position_V7(:,:)=data.eye_data(Ver,StartPoint:EndPoint,Select_Trial);%Vertical
    Position_H7=mean(Position_H7');
    Position_V7=mean(Position_V7');
    
    clear Select_Trial; Select_Trial=find(temp_stim_type==1 & temp_fp_rotate==3 & temp_rot_elevation==(i-1)*180+90);%FP_ROTATE=3 Trans only
    Position_H8(:,:)=data.eye_data(Hor,StartPoint:EndPoint, Select_Trial);%Horizontal
    Position_V8(:,:)=data.eye_data(Ver,StartPoint:EndPoint,Select_Trial);%Vertical
    Position_H8=mean(Position_H8');
    Position_V8=mean(Position_V8');
 


    figure(7);set(7,'Position', [5,5 980,650], 'name',FILE);orient landscape;% middle 2 sec
    subplot(4,2,i);plot(TimeLabel,Position_H5(:,:), 'r');
    title(['X/  Ro ',num2str((i-1)*180+90),'  (Roll)   Tra  ',num2str((i-1)*180+180),' deg / Tilt+Trans']);set(gca, 'xtick', [] );axis([1000 3000 13 18]);
    hold on;
    subplot(4,2,i+2);plot(TimeLabel,Position_H6(:,:), 'r');set(gca, 'xtick', [] );title(['Tilt-Trans']);axis([1000 3000 13 18]);hold on;
    subplot(4,2,i+4);plot(TimeLabel,Position_H7(:,:), 'r');set(gca, 'xtick', [] );title(['Tilt']);axis([1000 3000 13 18]);hold on;
    subplot(4,2,i+6);plot(TimeLabel,Position_H8(:,:), 'r');set(gca, 'xtick', [] );title(['Trans']);axis([1000 3000 13 18]);hold on;
    
    
    figure(8);set(8,'Position', [5,5 980,650], 'name',FILE);orient landscape;
    subplot(4,2,i);plot(TimeLabel,Position_V5(:,:), 'b');
    title(['Y/  Ro ',num2str((i-1)*180+90),'  (Roll)   Tra  ',num2str((i-1)*180+180),' deg / Tilt+Trans']);set(gca, 'xtick', [] );axis([1000 3000 13 18]);
    hold on;
    subplot(4,2,i+2);plot(TimeLabel,Position_V6(:,:), 'b');set(gca, 'xtick', [] );title(['Tilt-Trans']);axis([1000 3000 13 18]);hold on;
    subplot(4,2,i+4);plot(TimeLabel,Position_V7(:,:), 'b');set(gca, 'xtick', [] );title(['Tilt']);axis([1000 3000 13 18]);hold on;
    subplot(4,2,i+6);plot(TimeLabel,Position_V8(:,:), 'b');set(gca, 'xtick', [] );title(['Trans']);axis([1000 3000 13 18]);hold on;
    
    figure(9);set(9,'Position', [5,5 980,650], 'name',FILE);orient landscape;% whole 5 sec
    subplot(4,2,i);plot(TimeLabel,Position_H5(:,:), 'r');
    title(['X/  Ro ',num2str((i-1)*180+90),'  (Roll)   Tra  ',num2str((i-1)*180+180),' deg / Tilt+Trans']);set(gca, 'xtick', [] );
    hold on;
    subplot(4,2,i+2);plot(TimeLabel,Position_H6(:,:), 'r');set(gca, 'xtick', [] );title(['Tilt-Trans']);hold on;
    subplot(4,2,i+4);plot(TimeLabel,Position_H7(:,:), 'r');set(gca, 'xtick', [] );title(['Tilt']);hold on;
    subplot(4,2,i+6);plot(TimeLabel,Position_H8(:,:), 'r');set(gca, 'xtick', [] );title(['Trans']);hold on;
    
    
    figure(10);set(10,'Position', [5,5 980,650], 'name',FILE);orient landscape;
    subplot(4,2,i);plot(TimeLabel,Position_V5(:,:), 'b');
    title(['Y/  Ro ',num2str((i-1)*180+90),'  (Roll)   Tra  ',num2str((i-1)*180+180),' deg / Tilt+Trans']);set(gca, 'xtick', [] );
    hold on;
    subplot(4,2,i+2);plot(TimeLabel,Position_V6(:,:), 'b');set(gca, 'xtick', [] );title(['Tilt-Trans']);hold on;
    subplot(4,2,i+4);plot(TimeLabel,Position_V7(:,:), 'b');set(gca, 'xtick', [] );title(['Tilt']);hold on;
    subplot(4,2,i+6);plot(TimeLabel,Position_V8(:,:), 'b');set(gca, 'xtick', [] );title(['Trans']);hold on;
end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % output to text file
% % sprint_txt = ['%s'];
% % for i = 1 : 48
% %      sprint_txt = [sprint_txt, ' %1.2f'];    
% % end
% % %for j= 1:repeat
% % 	buff= sprintf(sprint_txt, FILE, azimuth_verg_LR(j,:) );
% %     %end
% % 	outfile = [BASE_PATH 'ProtocolSpecific\MOOG\3Dtuning\Eye_azimuthLR_Nor.dat'];
% % 	printflag = 0;
% % 	if (exist(outfile, 'file') == 0)    %file does not yet exist
% %         printflag = 1;
% % 	end
% % 	fid = fopen(outfile, 'a');
% % 	if (printflag)
% %         fprintf(fid, 'FILE\t');
% %         fprintf(fid, '\r\n');
% % 	end
% % for j=1:repeat
% % 	fprintf(fid, '%s', buff{j});
% %     fprintf(fid, '\r\n');
% % end
% % 	fclose(fid);

return;

