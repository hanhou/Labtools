function SimulatedData()

warning off MATLAB:divideByZero;
warning off MATLAB:singularMatrix;

clear;clc;close
timestep=0.1;
x = ([0 45 90 135 180 225 270 315 360]*pi/180)'; % for fitting
t = 0:timestep:2.5; 

clear global rawdata xdata tdata
global rawdata xdata tdata

for i = 1:length(t)
    xdata(:, i) = x;    
end
for i = 1:length(x)
    tdata(i,:) = t;
end
xtdata = [xdata; tdata];

%clear vect;vect=[13.4521  194.9922    1.2811    0    1.3468    0.2654   -0.1806    0.2415    180*pi/180];
clear vect;vect=[8.6315  126.0000    0.7374    4.9579    1.1834    0.1857   -0.4236    0.2000    0.0023];
clear vectVel;vectVel=[vect(1:7)    1    0];
% clear vectAcc1;vectAcc1=[vect(1:7)    0   vect(9)];
% clear vectAcc2;vectAcc2=[vect(1:3) 104*pi/180 vect(5:7) 0   0];
clear vectAcc1;vectAcc1=[8.6315  126.0000    0.8677    4.9579    1.1834    0.1857   -0.4236    0.2000    0.0023];
clear vectAcc2;vectAcc2=[8.6818  126.0000    0.8677    2.0529    1.1824    0.1832   -0.4006    0.2000    2.9051];

global model_use 
model_use=2;
spacetime_data = funccosnlin(vect,xtdata);rawdata=spacetime_data;
VelComp = vect(8)*funccosnlin(vectVel,xtdata);
AccComp1 = (1-vect(8))*funccosnlin(vectAcc1,xtdata);
AccComp2 = funccosnlin(vectAcc2,xtdata);

FigureIndex=2; 
figure(FigureIndex);set(FigureIndex,'Position', [50,50 1200,600], 'Name', 'CurveFitting');orient landscape; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot the rawdata
clear XAzi Ytime;[XAzi,Ytime] = meshgrid([0:45:360],[0:timestep:max(t)]);
figure(FigureIndex);axes('position',[0.05 0.77 0.17 0.15]);  
contourf(XAzi,Ytime,spacetime_data');colorbar;    
set(gca, 'xtick', [] );set(gca, 'XTickMode','manual'); 
set(gca, 'xtick',[0:90:360]);set(gca, 'xticklabel','0|90|180|270|360');
set(gca, 'ytick', [] ); set(gca, 'YTickMode','manual');
set(gca, 'ytick',[0:0.5:max(t)]);
set(gca, 'yticklabel','0|0.5|1|1.5|2|2.5');% set(gca, 'ydir','reverse');xlabel('Azimuth, X (deg)'); 
ylabel('Time (sec)');title('Raw Data');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot the velocity component
figure(FigureIndex);axes('position',[0.05 0.53 0.17 0.15]);
contourf(XAzi,Ytime,VelComp',10);colorbar;    
set(gca, 'xtick', [] );set(gca, 'XTickMode','manual'); 
set(gca, 'xtick',[0:90:315]);set(gca, 'xticklabel','0|90|180|270');
set(gca, 'ytick', [] ); set(gca, 'YTickMode','manual');
set(gca, 'ytick',[0:0.5:max(t)]);
set(gca, 'yticklabel','0|0.5|1|1.5|2|2.5');% set(gca, 'ydir','reverse');
%xlabel('Azimuth, X (deg)');  
ylabel('Time (sec)');title('Vel component  ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot the Acceleration component
figure(FigureIndex);axes('position',[0.05 0.31 0.17 0.15]);
contourf(XAzi,Ytime,AccComp1',10);colorbar;    
set(gca, 'xtick', [] );set(gca, 'XTickMode','manual'); 
set(gca, 'xtick',[0:90:315]);set(gca, 'xticklabel','0|90|180|270');
set(gca, 'ytick', [] ); set(gca, 'YTickMode','manual');
set(gca, 'ytick',[0:0.5:max(t)]);
set(gca, 'yticklabel','0|0.5|1|1.5|2|2.5');% set(gca, 'ydir','reverse');
% xlabel('Azimuth, X (deg)');  
ylabel('Time (sec)');title('Acc1 component  ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot the Acceleration component2
figure(FigureIndex);axes('position',[0.05 0.09 0.17 0.15]);
contourf(XAzi,Ytime,AccComp2',10);colorbar;    
set(gca, 'xtick', [] );set(gca, 'XTickMode','manual'); 
set(gca, 'xtick',[0:90:315]);set(gca, 'xticklabel','0|90|180|270');
set(gca, 'ytick', [] ); set(gca, 'YTickMode','manual');
set(gca, 'ytick',[0:0.5:max(t)]);
set(gca, 'yticklabel','0|0.5|1|1.5|2|2.5');% set(gca, 'ydir','reverse');
xlabel('Azimuth, X (deg)');  ylabel('Time (sec)');title('Acc2 component  ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out=[];
out = 'CurveFitting:    R0 |  Amp  |  n  | muAzi |  muT  | sigmaT |  DC2  | wVel  |ThetaAcc|  wAcc | ThetaPos | wPos'; 
out = strvcat(out, sprintf('------------------------------------------------------------------------------------------------------------'));    
clear OutputValue;OutputValue=vect;
OutputValue(4)=OutputValue(4)*180/pi;
OutputValue(9)=OutputValue(9)*180/pi;
out=strvcat(out, sprintf('Vel+Acc Model: %4.1f | %5.1f |%4.1f |%6.1f | %4.3f | %6.3f |%6.3f |%6.3f | %6.1f  ', OutputValue));

clear OutputValue;OutputValue=vectVel;
OutputValue(4)=OutputValue(4)*180/pi;
OutputValue(9)=OutputValue(9)*180/pi;
out=strvcat(out, sprintf('Vel Component: %4.1f | %5.1f |%4.1f |%6.1f | %4.3f | %6.3f |%6.3f |%6.3f | %6.1f  ', OutputValue));

clear OutputValue;OutputValue=vectAcc1;
OutputValue(4)=OutputValue(4)*180/pi;
OutputValue(9)=OutputValue(9)*180/pi;
out=strvcat(out, sprintf('AccComponent1: %4.1f | %5.1f |%4.1f |%6.1f | %4.3f | %6.3f |%6.3f |%6.3f | %6.1f  ', OutputValue));

clear OutputValue;OutputValue=vectAcc2;
OutputValue(4)=OutputValue(4)*180/pi;
OutputValue(9)=OutputValue(9)*180/pi;
out=strvcat(out, sprintf('AccComponent2: %4.1f | %5.1f |%4.1f |%6.1f | %4.3f | %6.3f |%6.3f |%6.3f | %6.1f  ', OutputValue));

figure(FigureIndex);axes('position',[0.26 0.69 0.25 0.25]); set(gca,'box','off','visible','off');
text(-0.08,1,out,'fontsize',8,'fontname','courier','horizontalalignment','left','verticalalignment','top');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% order in which the directions are plotted
plot_col = [1 1 1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 6 6 6 7 7 7 8 8 8];
plot_row = [5 4 3 2 1 4 3 2 4 3 2 4 3 2 4 3 2 4 3 2 4 3 2 4 3 2];

x_start = [0, 0];
x_stop =  [2, 2];
y_marker=[0,1.1*max(max(spacetime_data))];

xscale = [0.8 0.77 0.65 0.52 0.5 0.52 0.65 0.77];   
yscale = [0.35 0.5 0.65 0.5 0.35 0.2 0.05 0.2];

clear x_time;x_time=t;

UniAzi=[0:45:315];
for i=1:8
    i
    axes('position',[xscale(i) yscale(i) 0.16 0.10]);    
    plot(x_time,spacetime_data(i,:),'k','LineWidth',2);hold on;    
    plot(x_time,VelComp(i,:), 'r', 'LineWidth', 2);hold on;
    plot(x_time,AccComp1(i,:),'g','LineWidth',2);hold on;
    plot(x_time,AccComp2(i,:),'c','LineWidth',2);hold on;    
    
    xlim([0,2.5]);    ylim([0,1.2*max(max(spacetime_data))]);    
    set(gca, 'xtick', [] );set(gca, 'XTickMode','manual'); 
    set(gca, 'xtick',[0:0.5:2.5]);
    set( gca, 'xticklabel', '0|0.5|1|1.5|2|2.5');
    set( gca, 'yticklabel', ' ' );
    title({['azi = ' num2str(UniAzi(i))]}, 'FontSize', 8);
end  




