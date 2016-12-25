% analyze eyetrace data
clear all
% choose protocol
%
% analyze_protocol = [1];%Que
% analyze_protocol = [2];%Azrael
analyze_protocol = [3];%Zebulon
% analyze_protocol = [4];% Purusuit


if analyze_protocol == [1]
%     aa1 = dlmread('QueEye_Rot_Vet.dat','',1,1);%dim=size(aa1)  % load data
%     aa2 = dlmread('QueEye_Rot_Vis.dat','',1,1);%dim=size(aa2)  % load data
%     aa3 = dlmread('QueEye_Tra_Vet.dat','',1,1);%dim=size(aa3)  % load data
%     aa4 = dlmread('QueEye_Tra_Vis.dat','',1,1);%dim=size(aa4)  % load data
%     Labyrinthectomy 
    aa1 = dlmread('QueLaby_Rot_Vet.dat','',1,1);%dim=size(aa1)  % load data
    aa2 = dlmread('QueLaby_Rot_Vis.dat','',1,1);%dim=size(aa2)  % load data
    aa3 = dlmread('QueLaby_Tra_Vet.dat','',1,1);%dim=size(aa3)  % load data
    aa4 = dlmread('QueLaby_Tra_Vis.dat','',1,1);%dim=size(aa4)  % load data
   
elseif analyze_protocol == [2]
    aa1 = dlmread('AzraelEye_Rot_Vet.dat','',1,1);%dim=size(aa1)  % load data
    aa2 = dlmread('AzraelEye_Rot_Vis.dat','',1,1);%dim=size(aa2)  % load data
    aa3 = dlmread('AzraelEye_Tra_Vet.dat','',1,1);%dim=size(aa3)  % load data
    aa4 = dlmread('AzraelEye_Tra_Vis.dat','',1,1);%dim=size(aa4)  % load data
% %     aa = dlmread('Eye_rot_ves.dat','',1,1);  % load data
   
elseif analyze_protocol == [3]
%     aa3 = dlmread('ZebulonEye_Tra_vet.dat','',1,1);%dim=size(aa3)  % load data
%     aa4 = dlmread('ZebulonEye_Tra_vis.dat','',1,1);%dim=size(aa4)  % load data
%     aa1 = dlmread('ZebulonEye_Rot_vet.dat','',1,1);%dim=size(aa1)  % load data
%     aa2 = dlmread('ZebulonEye_Rot_vis.dat','',1,1);%dim=size(aa2)  % load data
%      Labyrinthectomy 
    aa1 = dlmread('ZebulonLaby_Rot_vet.dat','',1,1);%dim=size(aa1)  % load data
    aa2 = dlmread('ZebulonLaby_Rot_vis.dat','',1,1);%dim=size(aa2)  % load data
    aa3 = dlmread('ZebulonLaby_Tra_vet.dat','',1,1);%dim=size(aa3)  % load data
    aa4 = dlmread('ZebulonLaby_Tra_vis.dat','',1,1);%dim=size(aa4)  % load data
else
    aa = dlmread('Que_Pursuit_pursuit.dat','',1,1);  % load data 
%     aa = dlmread('Eye_rot_ves.dat','',1,1);  % load data
  
end

    title1 = ' Up    ';
    title2 = 'Down   ';
    title3 = 'Left   ';
    title4 = 'Right  ';
    
% mean for repeats cells...%Azuel no need
% aa3=aa3(:,1:22402);%Que 
% aa3=aa3(:,1:19202);%Zebulon or for Zebulon do not need

%%%%%%%%%%%%%%%%%%%%%% Select conditions !! %%%%%%%%%%%%%%%%
% aa=aa1;p_title = 'Rot Vest';%Rot_Vet
% aa=aa2;p_title = 'Rot Visu';%Rot_Vis
% aa=aa3;p_title = 'Tra Vest';%Tra_Vet
aa=aa4;p_title = 'Tra Visu';%Tra_Vis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








repeat = aa(:,1); % 1st column is repeatition, else are the raw eye trace
dim = size(aa)



for i = 1 : dim(1)  % How many cells?
   
        res_x_up(i,:) = aa(i, 2:401);
        res_y_up(i,:) = aa(i, 402:801);
        res_x_down(i,:)= aa(i, 802:1201);
        res_y_down(i,:)= aa(i, 1202:1601);
        res_x_left(i,:) = aa(i, 1602:2001);
        res_y_left(i,:) = aa(i, 2002:2401);
        res_x_right(i,:) = aa(i, 2402:2801);
        res_y_right(i,:) = aa(i, 2802:3201);
        
        
%         % Convert to velosity
%         vel_x_up(i,:) = diff(res_x_up(i,:))*1000/5;
% %         vel_x_up2{i}(j,:) = fderiv(res_x_up{i}(j,:),15,200);
%         vel_y_up(i,:) = diff(res_y_up(i,:))*200;
%         vel_x_down(i,:) = diff(res_x_down(i,:))*200;
%         vel_y_down(i,:) = diff(res_y_down(i,:))*200;
%         vel_x_left(i,:) = diff(res_x_left(i,:))*200;
%         vel_y_left(i,:) = diff(res_y_left(i,:))*200;
%         vel_x_right(i,:) = diff(res_x_right(i,:))*200;
%         vel_y_right(i,:) = diff(res_y_right(i,:))*200;
        
               % Convert to velosity 2
        vel_x_up2(i,:) = fderiv(res_x_up(i,:),15,200);
%         vel_x_up2{i}(j,:) = fderiv(res_x_up{i}(j,:),15,200);
        vel_y_up2(i,:) = fderiv(res_y_up(i,:),15,200);
        vel_x_down2(i,:) = fderiv(res_x_down(i,:),15,200);
        vel_y_down2(i,:) = fderiv(res_y_down(i,:),15,200);
        vel_x_left2(i,:) = fderiv(res_x_left(i,:),15,200);
        vel_y_left2(i,:) = fderiv(res_y_left(i,:),15,200);
        vel_x_right2(i,:) = fderiv(res_x_right(i,:),15,200);
        vel_y_right2(i,:) = fderiv(res_y_right(i,:),15,200);
end
        
  mean_p_up_101_299= mean((res_y_up(:,101:299))')';
  mean_p_down_101_299 = mean((res_y_down(:,101:299))')';
  mean_p_left_101_299 = mean((res_x_left(:,101:299))')';
  mean_p_right_101_299 = mean((res_x_right(:,101:299))')';

  mean_v2_up_101_299 = mean((vel_y_up2(:,101:299))')';
  mean_v2_down_101_299 = mean((vel_y_down2(:,101:299))')';
  mean_v2_left_101_299 = mean((vel_x_left2(:,101:299))')';
  mean_v2_right_101_299 = mean((vel_x_right2(:,101:299))')';

mp1sec=[mean_p_up_101_299 mean_p_down_101_299 mean_p_left_101_299 mean_p_right_101_299];
mv1sec=[mean_v2_up_101_299 mean_v2_down_101_299 mean_v2_left_101_299 mean_v2_right_101_299];

  mean_p_up_1_399 = mean((res_y_up(:,1:399))')';
  mean_p_down_1_399 = mean((res_y_down(:,1:399))')';
  mean_p_left_1_399 = mean((res_x_left(:,1:399))')';
  mean_p_right_1_399 = mean((res_x_right(:,1:399))')';
  
  mean_v2_up_1_399 = mean((vel_y_up2(:,1:399))')';
  mean_v2_down_1_399 = mean((vel_y_down2(:,1:399))')';
  mean_v2_left_1_399 = mean((vel_x_left2(:,1:399))')';
  mean_v2_right_1_399 = mean((vel_x_right2(:,1:399))')';        

mp2sec=[mean_p_up_1_399 mean_p_down_1_399 mean_p_left_1_399 mean_p_right_1_399];
mv2sec=[mean_v2_up_1_399 mean_v2_down_1_399 mean_v2_left_1_399 mean_v2_right_1_399];
  
  std_p_up_101_299= std((res_y_up(:,101:299))')';
  std_p_down_101_299 = std((res_y_down(:,101:299))')';
  std_p_left_101_299 = std((res_x_left(:,101:299))')';
  std_p_right_101_299 = std((res_x_right(:,101:299))')';
  
  
  std_v2_up_101_299 = std((vel_y_up2(:,101:299))')';
  std_v2_down_101_299 = std((vel_y_down2(:,101:299))')';
  std_v2_left_101_299 = std((vel_x_left2(:,101:299))')';
  std_v2_right_101_299 = std((vel_x_right2(:,101:299))')';

sp1sec=[std_p_up_101_299 std_p_down_101_299 std_p_left_101_299 std_p_right_101_299];
sv1sec=[std_v2_up_101_299 std_v2_down_101_299 std_v2_left_101_299 std_v2_right_101_299];

  std_p_up_1_399 = std((res_y_up(:,1:399))')';
  std_p_down_1_399 = std((res_y_down(:,1:399))')';
  std_p_left_1_399 = std((res_x_left(:,1:399))')';
  std_p_right_1_399 = std((res_x_right(:,1:399))')';
  
  std_v2_up_1_399 = std((vel_y_up2(:,1:399))')';
  std_v2_down_1_399 = std((vel_y_down2(:,1:399))')';
  std_v2_left_1_399 = std((vel_x_left2(:,1:399))')';
  std_v2_right_1_399 = std((vel_x_right2(:,1:399))')';  

sp2sec=[std_p_up_1_399 std_p_down_1_399 std_p_left_1_399 std_p_right_1_399];
sv2sec=[std_v2_up_1_399 std_v2_down_1_399 std_v2_left_1_399 std_v2_right_1_399];  
  
%%%%%%%%%%% %Save Data %%%%%%%%%%%%%%%%%5 
zero(1:dim(1))=[0];zero=zero';
ichi(1:dim(1))=[1];ichi=ichi';

onesec=[mp1sec zero sp1sec ichi mv1sec zero sv1sec];
twosec=[mp2sec zero sp2sec ichi mv2sec zero sv2sec];

csvwrite('onesec.dat',onesec);
csvwrite('twosec.dat',twosec);  
  
%%%%%%%%%%%%%%%%plot data%%%%%%%%%%%%%%%%%%%%%%
x=1:400;
i=[];

numfig=ceil(dim(1)/4); % 4 cells in 1 paper

for m=1:numfig  %% if paper limit... put the number

% for m=15:19    
    
figure(m+1);orient landscape;
% set(gca,'Position', [5,5 1000,680], 'Name', 'Envelope');    

 
i=(m-1)*4+1   %First 1 figure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1st

  subplot(4,4,1);

hl1=line(x, res_y_up(i,:),'Color','b','LineWidth',2);
ax1=gca;
set(ax1, 'XColor','k','YColor','b','Ylim',[-0.5 0.5])
set(ax1, 'XTickLabel', [])
%     set(gca, 'xtick',[1,100,200,300,400]);
%     set(gca, 'xticklabel','0|0.5|1|1.5|2'); 
set(get(ax1,'Ylabel'),'String','(deg)')
ax2=axes('Position', get(ax1,'Position'), 'Color', 'none', 'XColor', 'k', 'YColor', 'k', 'Ylim', [-3 3]);
set(ax2, 'YAxisLocation', 'right', 'ytick', [-2 0 2])
set(ax2, 'XAxisLocation', 'top', 'XTickLabel', [])
% set(get(ax2,'Ylabel'),'String','(deg/sec)')
hl2=line(x, vel_y_up2(i,:),'Color','k','Parent',ax2);
% title([num2str(mean_p_up_101_299(i)),num2str(std_p_up_101_299(i)),'/  ',title1,'/  ', p_title]);% '/  ',mean_v2_up_101_299(i),'+-',std_v2_up_101_299(i)]);
title([title1,'/  ', p_title, '      No.',num2str(i)]);    
    text (5,2,[num2str(mean_p_up_101_299(i)),' \pm ', num2str(std_p_up_101_299(i)), '  (deg)']);
    text (5,-2,[num2str(mean_v2_up_101_299(i)),' \pm ', num2str(std_v2_up_101_299(i)), '  (deg/s)']);
   
    
  subplot(4,4,5);
  
hl1=line(x, res_y_down(i,:),'Color','b','LineWidth',2);
ax1=gca;
set(ax1, 'XColor','k','YColor','b','Ylim',[-0.5 0.5])
set(ax1, 'XTickLabel', [])
%     set(gca, 'xtick',[1,100,200,300,400]);
%     set(gca, 'xticklabel','0|0.5|1|1.5|2'); 
set(get(ax1,'Ylabel'),'String','(deg)')
ax2=axes('Position', get(ax1,'Position'), 'Color', 'none', 'XColor', 'k', 'YColor', 'k', 'Ylim', [-3 3]);
set(ax2, 'YAxisLocation', 'right', 'ytick', [-2 0 2])
set(ax2, 'XAxisLocation', 'top', 'XTickLabel', [])
% set(get(ax2,'Ylabel'),'String','(deg/sec)')
hl2=line(x, vel_y_down2(i,:),'Color','k','Parent',ax2);
title([title2]);
    text (10,2,[num2str(mean_p_down_101_299(i)),' \pm ', num2str(std_p_down_101_299(i))]);
    text (20,-2,[num2str(mean_v2_down_101_299(i)),' \pm ', num2str(std_v2_down_101_299(i))]);
    
    
  subplot(4,4,9);
  
hl1=line(x, res_x_left(i,:),'Color','r','LineWidth',2);
ax1=gca;
set(ax1, 'XColor','k','YColor','b','Ylim',[-0.5 0.5])
set(ax1, 'XTickLabel', [])
%     set(gca, 'xtick',[1,100,200,300,400]);
%     set(gca, 'xticklabel','0|0.5|1|1.5|2'); 
set(get(ax1,'Ylabel'),'String','(deg)')
ax2=axes('Position', get(ax1,'Position'), 'Color', 'none', 'XColor', 'k', 'YColor', 'k', 'Ylim', [-3 3]);
set(ax2, 'YAxisLocation', 'right', 'ytick', [-2 0 2])
set(ax2, 'XAxisLocation', 'top', 'XTickLabel', [])
% set(get(ax2,'Ylabel'),'String','(deg/sec)')
hl2=line(x, vel_x_left2(i,:),'Color','k','Parent',ax2);
title([title3]);
    text (10,2,[num2str(mean_p_left_101_299(i)),' \pm ', num2str(std_p_left_101_299(i))]);
    text (20,-2,[num2str(mean_v2_left_101_299(i)),' \pm ', num2str(std_v2_left_101_299(i))]);
    
    
  subplot(4,4,13);
hl1=line(x, res_x_right(i,:),'Color','r','LineWidth',2);
ax1=gca;
set(ax1, 'XColor','k','YColor','b','Ylim',[-0.5 0.5])
% set(ax1, 'XTickLabel', [])
    set(gca, 'xtick',[1,100,200,300,400]);
    set(gca, 'xticklabel','0|0.5|1|1.5|2'); 
set(get(ax1,'Ylabel'),'String','(deg)')
set(get(ax1,'Xlabel'),'String','(sec)')
ax2=axes('Position', get(ax1,'Position'), 'Color', 'none', 'XColor', 'k', 'YColor', 'k', 'Ylim', [-3 3]);
set(ax2, 'YAxisLocation', 'right', 'ytick', [-2 0 2])
set(ax2, 'XAxisLocation', 'top', 'XTickLabel', [])
% set(get(ax2,'Ylabel'),'String','(deg/sec)')
hl2=line(x, vel_x_right2(i,:),'Color','k','Parent',ax2);
title([title4]);
    text (10,2,[num2str(mean_p_right_101_299(i)),' \pm ', num2str(std_p_right_101_299(i))]);
    text (20,-2,[num2str(mean_v2_right_101_299(i)),' \pm ', num2str(std_v2_right_101_299(i))]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2nd
i=i+1

   subplot(4,4,2);

hl1=line(x, res_y_up(i,:),'Color','b','LineWidth',2);
ax1=gca;
set(ax1, 'XColor','k','YColor','b','Ylim',[-0.5 0.5])
set(ax1, 'XTickLabel', [])
%     set(gca, 'xtick',[1,100,200,300,400]);
%     set(gca, 'xticklabel','0|0.5|1|1.5|2'); 
% set(get(ax1,'Ylabel'),'String','(deg)')
ax2=axes('Position', get(ax1,'Position'), 'Color', 'none', 'XColor', 'k', 'YColor', 'k', 'Ylim', [-3 3]);
set(ax2, 'YAxisLocation', 'right', 'ytick', [-2 0 2])
set(ax2, 'XAxisLocation', 'top', 'XTickLabel', [])
% set(get(ax2,'Ylabel'),'String','(deg/sec)')
hl2=line(x, vel_y_up2(i,:),'Color','k','Parent',ax2);
title([title1,'/  ', p_title, '      No.',num2str(i)]);    
    text (10,2,[num2str(mean_p_up_101_299(i)),' \pm ', num2str(std_p_up_101_299(i))]);
    text (20,-2,[num2str(mean_v2_up_101_299(i)),' \pm ', num2str(std_v2_up_101_299(i))]);
   
    
  subplot(4,4,6);
  
hl1=line(x, res_y_down(i,:),'Color','b','LineWidth',2);
ax1=gca;
set(ax1, 'XColor','k','YColor','b','Ylim',[-0.5 0.5])
set(ax1, 'XTickLabel', [])
%     set(gca, 'xtick',[1,100,200,300,400]);
%     set(gca, 'xticklabel','0|0.5|1|1.5|2'); 
% set(get(ax1,'Ylabel'),'String','(deg)')
ax2=axes('Position', get(ax1,'Position'), 'Color', 'none', 'XColor', 'k', 'YColor', 'k', 'Ylim', [-3 3]);
set(ax2, 'YAxisLocation', 'right', 'ytick', [-2 0 2])
set(ax2, 'XAxisLocation', 'top', 'XTickLabel', [])
% set(get(ax2,'Ylabel'),'String','(deg/sec)')
hl2=line(x, vel_y_down2(i,:),'Color','k','Parent',ax2);
title([title2]);
    text (10,2,[num2str(mean_p_down_101_299(i)),' \pm ', num2str(std_p_down_101_299(i))]);
    text (20,-2,[num2str(mean_v2_down_101_299(i)),' \pm ', num2str(std_v2_down_101_299(i))]);
 
    
    
  subplot(4,4,10);
  
hl1=line(x, res_x_left(i,:),'Color','r','LineWidth',2);
ax1=gca;
set(ax1, 'XColor','k','YColor','b','Ylim',[-0.5 0.5])
set(ax1, 'XTickLabel', [])
%     set(gca, 'xtick',[1,100,200,300,400]);
%     set(gca, 'xticklabel','0|0.5|1|1.5|2'); 
% set(get(ax1,'Ylabel'),'String','(deg)')
ax2=axes('Position', get(ax1,'Position'), 'Color', 'none', 'XColor', 'k', 'YColor', 'k', 'Ylim', [-3 3]);
set(ax2, 'YAxisLocation', 'right', 'ytick', [-2 0 2])
set(ax2, 'XAxisLocation', 'top', 'XTickLabel', [])
% set(get(ax2,'Ylabel'),'String','(deg/sec)')
hl2=line(x, vel_x_left2(i,:),'Color','k','Parent',ax2);
title([title3]);
    text (10,2,[num2str(mean_p_left_101_299(i)),' \pm ', num2str(std_p_left_101_299(i))]);
    text (20,-2,[num2str(mean_v2_left_101_299(i)),' \pm ', num2str(std_v2_left_101_299(i))]);
  
    
    
  subplot(4,4,14);
  
hl1=line(x, res_x_right(i,:),'Color','r','LineWidth',2);
ax1=gca;
set(ax1, 'XColor','k','YColor','b','Ylim',[-0.5 0.5])
% set(ax1, 'XTickLabel', [])
    set(gca, 'xtick',[1,100,200,300,400]);
    set(gca, 'xticklabel','0|0.5|1|1.5|2'); 
% set(get(ax1,'Ylabel'),'String','(deg)')
ax2=axes('Position', get(ax1,'Position'), 'Color', 'none', 'XColor', 'k', 'YColor', 'k', 'Ylim', [-3 3]);
set(ax2, 'YAxisLocation', 'right', 'ytick', [-2 0 2])
set(ax2, 'XAxisLocation', 'top', 'XTickLabel', [])
% set(get(ax2,'Ylabel'),'String','(deg/sec)')
hl2=line(x, vel_x_right2(i,:),'Color','k','Parent',ax2);
title([title4]);
    text (10,2,[num2str(mean_p_right_101_299(i)),' \pm ', num2str(std_p_right_101_299(i))]);
    text (20,-2,[num2str(mean_v2_right_101_299(i)),' \pm ', num2str(std_v2_right_101_299(i))]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 3rd
i=i+1

  subplot(4,4,3);

hl1=line(x, res_y_up(i,:),'Color','b','LineWidth',2);
ax1=gca;
set(ax1, 'XColor','k','YColor','b','Ylim',[-0.5 0.5])
set(ax1, 'XTickLabel', [])
%     set(gca, 'xtick',[1,100,200,300,400]);
%     set(gca, 'xticklabel','0|0.5|1|1.5|2'); 
% set(get(ax1,'Ylabel'),'String','(deg)')
ax2=axes('Position', get(ax1,'Position'), 'Color', 'none', 'XColor', 'k', 'YColor', 'k', 'Ylim', [-3 3]);
set(ax2, 'YAxisLocation', 'right', 'ytick', [-2 0 2])
set(ax2, 'XAxisLocation', 'top', 'XTickLabel', [])
% set(get(ax2,'Ylabel'),'String','(deg/sec)')
hl2=line(x, vel_y_up2(i,:),'Color','k','Parent',ax2);
title([title1,'/  ', p_title, '      No.',num2str(i)]);    
    text (10,2,[num2str(mean_p_up_101_299(i)),' \pm ', num2str(std_p_up_101_299(i))]);
    text (20,-2,[num2str(mean_v2_up_101_299(i)),' \pm ', num2str(std_v2_up_101_299(i))]);
   
    
  subplot(4,4,7);
  
hl1=line(x, res_y_down(i,:),'Color','b','LineWidth',2);
ax1=gca;
set(ax1, 'XColor','k','YColor','b','Ylim',[-0.5 0.5])
set(ax1, 'XTickLabel', [])
%     set(gca, 'xtick',[1,100,200,300,400]);
%     set(gca, 'xticklabel','0|0.5|1|1.5|2'); 
% set(get(ax1,'Ylabel'),'String','(deg)')
ax2=axes('Position', get(ax1,'Position'), 'Color', 'none', 'XColor', 'k', 'YColor', 'k', 'Ylim', [-3 3]);
set(ax2, 'YAxisLocation', 'right', 'ytick', [-2 0 2])
set(ax2, 'XAxisLocation', 'top', 'XTickLabel', [])
% set(get(ax2,'Ylabel'),'String','(deg/sec)')
hl2=line(x, vel_y_down2(i,:),'Color','k','Parent',ax2);
title([title2]);
    text (10,2,[num2str(mean_p_down_101_299(i)),' \pm ', num2str(std_p_down_101_299(i))]);
    text (20,-2,[num2str(mean_v2_down_101_299(i)),' \pm ', num2str(std_v2_down_101_299(i))]);
 
    
    
  subplot(4,4,11);
  
hl1=line(x, res_x_left(i,:),'Color','r','LineWidth',2);
ax1=gca;
set(ax1, 'XColor','k','YColor','b','Ylim',[-0.5 0.5])
set(ax1, 'XTickLabel', [])
%     set(gca, 'xtick',[1,100,200,300,400]);
%     set(gca, 'xticklabel','0|0.5|1|1.5|2'); 
% set(get(ax1,'Ylabel'),'String','(deg)')
ax2=axes('Position', get(ax1,'Position'), 'Color', 'none', 'XColor', 'k', 'YColor', 'k', 'Ylim', [-3 3]);
set(ax2, 'YAxisLocation', 'right', 'ytick', [-2 0 2])
set(ax2, 'XAxisLocation', 'top', 'XTickLabel', [])
% set(get(ax2,'Ylabel'),'String','(deg/sec)')
hl2=line(x, vel_x_left2(i,:),'Color','k','Parent',ax2);
title([title3]);
    text (10,2,[num2str(mean_p_left_101_299(i)),' \pm ', num2str(std_p_left_101_299(i))]);
    text (20,-2,[num2str(mean_v2_left_101_299(i)),' \pm ', num2str(std_v2_left_101_299(i))]);
  
    
    
  subplot(4,4,15);
  
hl1=line(x, res_x_right(i,:),'Color','r','LineWidth',2);
ax1=gca;
set(ax1, 'XColor','k','YColor','b','Ylim',[-0.5 0.5])
% set(ax1, 'XTickLabel', [])
    set(gca, 'xtick',[1,100,200,300,400]);
    set(gca, 'xticklabel','0|0.5|1|1.5|2'); 
% set(get(ax1,'Ylabel'),'String','(deg)')
ax2=axes('Position', get(ax1,'Position'), 'Color', 'none', 'XColor', 'k', 'YColor', 'k', 'Ylim', [-3 3]);
set(ax2, 'YAxisLocation', 'right', 'ytick', [-2 0 2])
set(ax2, 'XAxisLocation', 'top', 'XTickLabel', [])
% set(get(ax2,'Ylabel'),'String','(deg/sec)')
hl2=line(x, vel_x_right2(i,:),'Color','k','Parent',ax2);
title([title4]);
    text (10,2,[num2str(mean_p_right_101_299(i)),' \pm ', num2str(std_p_right_101_299(i))]);
    text (20,-2,[num2str(mean_v2_right_101_299(i)),' \pm ', num2str(std_v2_right_101_299(i))]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 4th
i=i+1

   subplot(4,4,4);

hl1=line(x, res_y_up(i,:),'Color','b','LineWidth',2);
ax1=gca;
set(ax1, 'XColor','k','YColor','b','Ylim',[-0.5 0.5])
set(ax1, 'XTickLabel', [])
%     set(gca, 'xtick',[1,100,200,300,400]);
%     set(gca, 'xticklabel','0|0.5|1|1.5|2'); 
% set(get(ax1,'Ylabel'),'String','(deg)')
ax2=axes('Position', get(ax1,'Position'), 'Color', 'none', 'XColor', 'k', 'YColor', 'k', 'Ylim', [-3 3]);
set(ax2, 'YAxisLocation', 'right', 'ytick', [-2 0 2])
set(ax2, 'XAxisLocation', 'top', 'XTickLabel', [])
set(get(ax2,'Ylabel'),'String','(deg/sec)')
hl2=line(x, vel_y_up2(i,:),'Color','k','Parent',ax2);
title([title1,'/  ', p_title, '      No.',num2str(i)]);    
    text (10,2,[num2str(mean_p_up_101_299(i)),' \pm ', num2str(std_p_up_101_299(i))]);
    text (20,-2,[num2str(mean_v2_up_101_299(i)),' \pm ', num2str(std_v2_up_101_299(i))]);
   
    
  subplot(4,4,8);
  
hl1=line(x, res_y_down(i,:),'Color','b','LineWidth',2);
ax1=gca;
set(ax1, 'XColor','k','YColor','b','Ylim',[-0.5 0.5])
set(ax1, 'XTickLabel', [])
%     set(gca, 'xtick',[1,100,200,300,400]);
%     set(gca, 'xticklabel','0|0.5|1|1.5|2'); 
% set(get(ax1,'Ylabel'),'String','(deg)')
ax2=axes('Position', get(ax1,'Position'), 'Color', 'none', 'XColor', 'k', 'YColor', 'k', 'Ylim', [-3 3]);
set(ax2, 'YAxisLocation', 'right', 'ytick', [-2 0 2])
set(ax2, 'XAxisLocation', 'top', 'XTickLabel', [])
set(get(ax2,'Ylabel'),'String','(deg/sec)')
hl2=line(x, vel_y_down2(i,:),'Color','k','Parent',ax2);
title([title2]);
    text (10,2,[num2str(mean_p_down_101_299(i)),' \pm ', num2str(std_p_down_101_299(i))]);
    text (20,-2,[num2str(mean_v2_down_101_299(i)),' \pm ', num2str(std_v2_down_101_299(i))]);
 
    
    
  subplot(4,4,12);
  
hl1=line(x, res_x_left(i,:),'Color','r','LineWidth',2);
ax1=gca;
set(ax1, 'XColor','k','YColor','b','Ylim',[-0.5 0.5])
set(ax1, 'XTickLabel', [])
%     set(gca, 'xtick',[1,100,200,300,400]);
%     set(gca, 'xticklabel','0|0.5|1|1.5|2'); 
% set(get(ax1,'Ylabel'),'String','(deg)')
ax2=axes('Position', get(ax1,'Position'), 'Color', 'none', 'XColor', 'k', 'YColor', 'k', 'Ylim', [-3 3]);
set(ax2, 'YAxisLocation', 'right', 'ytick', [-2 0 2])
set(ax2, 'XAxisLocation', 'top', 'XTickLabel', [])
set(get(ax2,'Ylabel'),'String','(deg/sec)')
hl2=line(x, vel_x_left2(i,:),'Color','k','Parent',ax2);
title([title3]);
    text (10,2,[num2str(mean_p_left_101_299(i)),' \pm ', num2str(std_p_left_101_299(i))]);
    text (20,-2,[num2str(mean_v2_left_101_299(i)),' \pm ', num2str(std_v2_left_101_299(i))]);
  
    
    
  subplot(4,4,16);
  
hl1=line(x, res_x_right(i,:),'Color','r','LineWidth',2);
ax1=gca;
set(ax1, 'XColor','k','YColor','b','Ylim',[-0.5 0.5])
% set(ax1, 'XTickLabel', [])
    set(gca, 'xtick',[1,100,200,300,400]);
    set(gca, 'xticklabel','0|0.5|1|1.5|2'); 
% set(get(ax1,'Ylabel'),'String','(deg)')
ax2=axes('Position', get(ax1,'Position'), 'Color', 'none', 'XColor', 'k', 'YColor', 'k', 'Ylim', [-3 3]);
set(ax2, 'YAxisLocation', 'right', 'ytick', [-2 0 2])
set(ax2, 'XAxisLocation', 'top', 'XTickLabel', [])
set(get(ax2,'Ylabel'),'String','(deg/sec)')
hl2=line(x, vel_x_right2(i,:),'Color','k','Parent',ax2);
title([title4]);
    text (10,2,[num2str(mean_p_right_101_299(i)),' \pm ', num2str(std_p_right_101_299(i))]);
    text (20,-2,[num2str(mean_v2_right_101_299(i)),' \pm ', num2str(std_v2_right_101_299(i))]);


i=[];%clear i, calculated by m

end       



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% % middle 1 sec 101-299
% % all 2 sec 1 - 399
% %
% 
% 
% % 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5        output files  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% p_up_1=mean(p_up_101_299);
% p_down_1=mean(p_down_101_299);
% p_left_1=mean(p_left_101_299);
% p_right_1=mean(p_right_101_299);
% 
% v1_up_1=mean(v1_up_101_299);
% v1_down_1=mean(v1_down_101_299);
% v1_left_1=mean(v1_left_101_299);
% v1_right_1=mean(v1_right_101_299);
% 
% v2_up_1=mean(v2_up_101_299);
% v2_down_1=mean(v2_down_101_299);
% v2_left_1=mean(v2_left_101_299);
% v2_right_1=mean(v2_right_101_299);
% 
%     space=[0];
% 
% p_up_2=mean(p_up_1_399);
% p_down_2=mean(p_down_1_399);
% p_left_2=mean(p_left_1_399);
% p_right_2=mean(p_right_1_399);
% 
% v1_up_2=mean(v1_up_1_399);
% v1_down_2=mean(v1_down_1_399);
% v1_left_2=mean(v1_left_1_399);
% v1_right_2=mean(v1_right_1_399);
% 
% v2_up_2=mean(v2_up_1_399);
% v2_down_2=mean(v2_down_1_399);
% v2_left_2=mean(v2_left_1_399);
% v2_right_2=mean(v2_right_1_399);   
%     
% middle_mean=[p_up_1 p_down_1 p_left_1 p_right_1 space v1_up_1 v1_down_1 v1_left_1 v1_right_1 space v2_up_1 v2_down_1 v2_left_1 v2_right_1];
% whole_mean=[p_up_2 p_down_2 p_left_2 p_right_2 space v1_up_2 v1_down_2 v1_left_2 v1_right_2 space v2_up_2 v2_down_2 v2_left_2 v2_right_2];
% 
% %%%%%%%%%%%%%%%%%%%% STD %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sp_up_1=std(p_up_101_299);
% sp_down_1=std(p_down_101_299);
% sp_left_1=std(p_left_101_299);
% sp_right_1=std(p_right_101_299);
% 
% sv1_up_1=std(v1_up_101_299);
% sv1_down_1=std(v1_down_101_299);
% sv1_left_1=std(v1_left_101_299);
% sv1_right_1=std(v1_right_101_299);
% 
% sv2_up_1=std(v2_up_101_299);
% sv2_down_1=std(v2_down_101_299);
% sv2_left_1=std(v2_left_101_299);
% sv2_right_1=std(v2_right_101_299);
% 
% %     space=[0];
% 
% sp_up_2=std(p_up_1_399);
% sp_down_2=std(p_down_1_399);
% sp_left_2=std(p_left_1_399);
% sp_right_2=std(p_right_1_399);
% 
% sv1_up_2=std(v1_up_1_399);
% sv1_down_2=std(v1_down_1_399);
% sv1_left_2=std(v1_left_1_399);
% sv1_right_2=std(v1_right_1_399);
% 
% sv2_up_2=std(v2_up_1_399);
% sv2_down_2=std(v2_down_1_399);
% sv2_left_2=std(v2_left_1_399);
% sv2_right_2=std(v2_right_1_399);   
%     
% middle_std=[sp_up_1 sp_down_1 sp_left_1 sp_right_1 space sv1_up_1 sv1_down_1 sv1_left_1 sv1_right_1 space sv2_up_1 sv2_down_1 sv2_left_1 sv2_right_1];
% whole_std=[sp_up_2 sp_down_2 sp_left_2 sp_right_2 space sv1_up_2 sv1_down_2 sv1_left_2 sv1_right_2 space sv2_up_2 sv2_down_2 sv2_left_2 sv2_right_2];
% 
% 
%     summary=[middle_mean;middle_std;whole_mean;whole_std]
%     
%     
%     
%     
% %     csvwrite('summary_101_299.dat',summary_101_299);
%     
