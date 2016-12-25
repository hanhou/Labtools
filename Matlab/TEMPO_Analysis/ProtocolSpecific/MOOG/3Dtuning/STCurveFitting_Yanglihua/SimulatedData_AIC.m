warning off MATLAB:divideByZero;
warning off MATLAB:singularMatrix;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;clc
timestep=0.1;
x = ([0 45 90 135 180 225 270 315 360]*pi/180)'; % for fitting
t = 0:timestep:2.5; 

for i = 1:length(t)
    xdata(:, i) = x;    
end
for i = 1:length(x)
    tdata(i,:) = t;
end
xtdata = [xdata; tdata];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%VelOnly model
% vect=[13.4 106.7 1.2 191.4/180*pi 1.349 0.273 0.664 0.999 170.0/180*pi];
% vect=[13.4 106.7 1.2 191.4/180*pi 1.349 0.273 0.664 0.253 170.0/180*pi];
vect=[13.4 114.7 1.0 176.0/180*pi 1.302 0.261 0.627 0.01 194.2/180*pi];

%Position model
% vect=[12.14 133.1 1.2 183.8/180*pi 1.290 0.254 0.566 0.708 86.3/180*pi 238.0*pi/180 0.851];
% vect=[12.14 50.1 1.2 183.8/180*pi 1.290 0.254 0.566 0.708 86.3/180*pi 238.0*pi/180 0.851];
% vect=[12.14 50.1 1.2 183.8/180*pi 1.290 0.254 0.566 0.200 86.3/180*pi 238.0*pi/180 0.851];
% vect=[13.1 127.0 1.1 191.8/180*pi 1.289 0.259 0.608 0.225 170.4/180*pi 219.2*pi/180 0.018];
% vect=[13.1 127.0 1.1 191.8/180*pi 1.289 0.259 0.608 0.1 170.4/180*pi 0.01 219.2*pi/180 0.9];
% vect=[13.1 127.0 1.1 191.8/180*pi 1.289 0.259 0.608 0.1 170.4/180*pi 0.01 219.2*pi/180 0.9];
% vect=[12.1 133.049 1.18 183.8/180*pi 1.290 0.254 0.566 0.308 86.3/180*pi 0.014 238.0*pi/180 0.6780];
% vect=[12.1 133.049 1.18 183.8/180*pi 1.290 0.254 0.566 0.308 86.3/180*pi 0.014 238.0*pi/180 0.6780];

global model_use
model_use=1;
spacetime_data = funccosnlin(vect,xtdata);
desired_resp=reshape(spacetime_data,1,size(spacetime_data,1)*size(spacetime_data,2));
Nreps=5;
% figure;contourf(spacetime_data');colorbar;

%loop through a bunch of randomized simulations
%now generate one set of fake data by Poisson sampling
for i=1:length(desired_resp)%length(disp)
    i;
    rand_resps = poissrnd(desired_resp(i)*ones(Nreps));
    for j=1:Nreps
        plot_y(j,i) = rand_resps(j);
    end
end
means=mean(plot_y);means=reshape(means,size(spacetime_data,1),size(spacetime_data,2));
FigureIndex=2; 
figure(FigureIndex);set(FigureIndex,'Position', [50,100 1200,800], 'Name', 'CurveFitting');orient landscape; 
clear XAzi Ytime;[XAzi,Ytime] = meshgrid([0:45:360],[0:timestep:max(t)]);
figure(FigureIndex);axes('position',[0.05 0.77 0.17 0.15]);   
contourf(XAzi,Ytime,means');colorbar;    
set(gca, 'xtick', [] );set(gca, 'XTickMode','manual'); 
set(gca, 'xtick',[0:90:360]);set(gca, 'xticklabel','0|90|180|270|360');
set(gca, 'ytick', [] ); set(gca, 'YTickMode','manual');
set(gca, 'ytick',[0:0.5:max(t)]);
set(gca, 'yticklabel','0|0.5|1|1.5|2|2.5');% set(gca, 'ydir','reverse');
xlabel('Azimuth, X (deg)');  ylabel('Time (sec)');
title('Raw data');

%------------Velocity FITS----------------- 
global model_use
model_use=1;
[spacefit_Vel,vect_Vel,r_squared_Vel,CI_Vel] = MSF_Vel_fit(means,timestep,1);    
figure(FigureIndex);axes('position',[0.05 0.53 0.17 0.15]);
contourf(XAzi,Ytime,spacefit_Vel',10);colorbar;    
set(gca, 'xtick', [] );set(gca, 'XTickMode','manual'); 
set(gca, 'xtick',[0:90:360]);set(gca, 'xticklabel','0|90|180|270|360');
set(gca, 'ytick', [] ); set(gca, 'YTickMode','manual');
set(gca, 'ytick',[0:0.5:max(t)]);
set(gca, 'yticklabel','0|0.5|1|1.5|2|2.5');% set(gca, 'ydir','reverse');%xlabel('Azimuth, X (deg)');  
ylabel('Time (sec)');title('model: Vel  ');

error_surf = means - spacefit_Vel;
err_Vel = cosnlin_err(vect_Vel);
figure(FigureIndex);axes('position',[0.26 0.53 0.17 0.15]);%subplot('position', [0.1 0.7 0.22 0.22]);
contourf(XAzi,Ytime,error_surf'); colorbar;    
set(gca, 'xtick', [] );set(gca, 'XTickMode','manual'); 
set(gca, 'xtick',[0:90:360]);set(gca, 'xticklabel','0|90|180|270|360');
set(gca, 'ytick', [] ); set(gca, 'YTickMode','manual');
set(gca, 'ytick',[0:0.5:max(t)]);
set(gca, 'yticklabel','0|0.5|1|1.5|2|2.5');% set(gca, 'ydir','reverse');
title({[ 'Err: ' num2str(err_Vel, '%0.2f') ]},  'FontSize', 10);

%Do the chi-square goodness of fit test
[chi2_Vel, chi2P_Vel] = Chi2_Test_3D(xtdata, plot_y, 'funccosnlin', vect_Vel, length(vect_Vel));

%------------Velocity+Acceleration FITS-----------------    
model_use=2;
[spacefit_VelAcc,vect_VelAcc,r_squared_VelAcc,CI_VelAcc] = MSF_VelAcc_fit(means,vect_Vel,timestep,1);        
figure(FigureIndex);axes('position',[0.05 0.31 0.17 0.15]);
contourf(XAzi,Ytime,spacefit_VelAcc');colorbar;    
set(gca, 'xtick', [] );set(gca, 'XTickMode','manual'); 
set(gca, 'xtick',[0:90:360]);set(gca, 'xticklabel','0|90|180|270|360');
set(gca, 'ytick', [] ); set(gca, 'YTickMode','manual');
set(gca, 'ytick',[0:0.5:max(t)]);
set(gca, 'yticklabel','0|0.5|1|1.5|2|2.5');% set(gca, 'ydir','reverse');
%xlabel('Azimuth, X (deg)');  
ylabel('Time (sec)');title('model: Vel + Acc ');

error_surf = means - spacefit_VelAcc;%current_err = cosnlin_err(vect_VelAcc);
err_VelAcc = cosnlin_err(vect_VelAcc);
figure(FigureIndex);axes('position',[0.26 0.31 0.17 0.15]);
contourf(XAzi,Ytime,error_surf'); colorbar;   
set(gca, 'xtick', [] );set(gca, 'XTickMode','manual'); 
set(gca, 'xtick',[0:90:360]);set(gca, 'xticklabel','0|90|180|270|360');
set(gca, 'ytick', [] ); set(gca, 'YTickMode','manual');
set(gca, 'ytick',[0:0.5:max(t)]);
set(gca, 'yticklabel','0|0.5|1|1.5|2|2.5');% set(gca, 'ydir','reverse');
title({[ 'Err: ' num2str(err_VelAcc, '%0.2f') ]},  'FontSize', 10);% axis off;

%Do the chi-square goodness of fit test
[chi2_VelAcc, chi2P_VelAcc] = Chi2_Test_3D(xtdata, plot_y, 'funccosnlin', vect_VelAcc, length(vect_VelAcc));

%------------Velocity+Acc+Position FITS-----------------    
model_use=3;
[spacefit_VelAccPos,vect_VelAccPos,r_squared_VelAccPos,CI_VelAccPos] = MSF_VelAccPos_fit(means,vect_VelAcc,timestep,1);        
figure(FigureIndex);axes('position',[0.05 0.09 0.17 0.15]);
contourf(XAzi,Ytime,spacefit_VelAccPos');colorbar;    
set(gca, 'xtick', [] );set(gca, 'XTickMode','manual'); 
set(gca, 'xtick',[0:90:360]);set(gca, 'xticklabel','0|90|180|270|360');
set(gca, 'ytick', [] ); set(gca, 'YTickMode','manual');
set(gca, 'ytick',[0:0.5:max(t)]);
set(gca, 'yticklabel','0|0.5|1|1.5|2|2.5');% set(gca, 'ydir','reverse');
%xlabel('Azimuth, X (deg)');  
ylabel('Time (sec)');title('model: Vel + Acc + Pos ');

error_surf = means - spacefit_VelAccPos;%current_err = cosnlin_err(vect_VelAcc);
err_VelAccPos = cosnlin_err(vect_VelAccPos);
figure(FigureIndex);axes('position',[0.26 0.09 0.17 0.15]);
contourf(XAzi,Ytime,error_surf'); colorbar;   
set(gca, 'xtick', [] );set(gca, 'XTickMode','manual'); 
set(gca, 'xtick',[0:90:360]);set(gca, 'xticklabel','0|90|180|270|360');
set(gca, 'ytick', [] ); set(gca, 'YTickMode','manual');
set(gca, 'ytick',[0:0.5:max(t)]);
set(gca, 'yticklabel','0|0.5|1|1.5|2|2.5');% set(gca, 'ydir','reverse');
title({[ 'Err: ' num2str(err_VelAccPos, '%0.2f') ]},  'FontSize', 10);% axis off;

%Do the chi-square goodness of fit test
[chi2_VelAccPos, chi2P_VelAccPos] = Chi2_Test_3D(xtdata, plot_y, 'funccosnlin', vect_VelAccPos, length(vect_VelAccPos));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out=[];
out = 'CurveFitting:    R0 |  Amp  |  n  | muAzi |  muT  | sigmaT |  DC2  | wVel  |ThetaAcc|  wAcc | ThetaPos | wPos'; 
out = strvcat(out, sprintf('------------------------------------------------------------------------------------------------------------'));    
clear OutputValue;OutputValue=vect_Vel;
OutputValue(4)=OutputValue(4)*180/pi;
out=strvcat(out, sprintf('Vel Model:     %4.1f | %5.1f |%4.1f |%6.1f | %4.3f | %6.3f |%6.3f', OutputValue));

clear OutputValue;OutputValue=vect_VelAcc;
OutputValue(4)=OutputValue(4)*180/pi;
OutputValue(9)=OutputValue(9)*180/pi;
out=strvcat(out, sprintf('Vel+Acc Model: %4.1f | %5.1f |%4.1f |%6.1f | %4.3f | %6.3f |%6.3f |%6.3f | %6.1f  ', OutputValue));

clear OutputValue; OutputValue=[vect_VelAccPos(1:9) 1-vect_VelAccPos(8) vect_VelAccPos(10:11)];
OutputValue(4)=OutputValue(4)*180/pi;
OutputValue(9)=OutputValue(9)*180/pi;
OutputValue(11)=OutputValue(11)*180/pi;
out=strvcat(out, sprintf('Vel+Acc+Pos  : %4.1f | %5.1f |%4.1f |%6.1f | %4.3f | %6.3f |%6.3f |%6.3f | %6.1f |%6.3f | %6.1f   |%6.3f', OutputValue)); 
figure(FigureIndex);axes('position',[0.26 0.74 0.25 0.25]); set(gca,'box','off','visible','off');
text(-0.08,1,out,'fontsize',8,'fontname','courier','horizontalalignment','left','verticalalignment','top');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do sequential F-test
Ftest_1vs2=[(err_Vel-err_VelAcc)/(length(vect_VelAcc)-length(vect_Vel))]/[err_VelAcc/(size(spacetime_data,1)*size(spacetime_data,2)-length(vect_Vel))];
p_1vs2=1-fcdf(Ftest_1vs2,length(vect_VelAcc)-length(vect_Vel),size(spacetime_data,1)*size(spacetime_data,2)-length(vect_Vel));

Ftest_2vs3=[(err_VelAcc-err_VelAccPos)/(length(vect_VelAccPos)-length(vect_VelAcc))]/[err_VelAccPos/(size(spacetime_data,1)*size(spacetime_data,2)-length(vect_VelAcc))];
p_2vs3=1-fcdf(Ftest_2vs3,length(vect_VelAccPos)-length(vect_VelAcc),size(spacetime_data,1)*size(spacetime_data,2)-length(vect_VelAcc));

% do AIC test
AIC_1vs2=size(spacetime_data,1)*size(spacetime_data,2)*log(err_VelAcc/err_Vel)+2*(length(vect_VelAcc)-length(vect_Vel));
AIC_2vs3=size(spacetime_data,1)*size(spacetime_data,2)*log(err_VelAccPos/err_VelAcc)+2*(length(vect_VelAccPos)-length(vect_VelAcc));

out=[];
% out = 'Goodness of fitting:  r2 (Vel) | r2 (VelAcc)| r2 (VelAccPos) | Ftest(1vs2) | p(1vs2)  | Ftest(2vs3) | p(2vs3)'; 
out = 'Goodness of fit: r2(Vel) | r2(V+A)| r2(V+A+P) | F(1vs2) | p(1vs2)  | F(2vs3) | p(2vs3) | AIC(1vs2) | AIC(2vs3)'; 
out = strvcat(out, sprintf('--------------------------------------------------------------------------------------------------------'));    
clear OutputValue;OutputValue=[r_squared_Vel r_squared_VelAcc r_squared_VelAccPos Ftest_1vs2 p_1vs2 Ftest_2vs3 p_2vs3 AIC_1vs2 AIC_2vs3];
out=strvcat(out, sprintf('                  %6.3f | %6.3f |  %6.3f   | %6.3f  |  %6.3f  |  %6.3f |  %6.3f |    %6.3f |  %6.3f ', OutputValue));
figure(FigureIndex);axes('position',[0.26 0.65 0.25 0.25]); set(gca,'box','off','visible','off');
text(-0.08,1,out,'fontsize',8,'fontname','courier','horizontalalignment','left','verticalalignment','top');

out=[];
out = 'Goodness of fit: Chi2(Vel) | p(Vel) | Chi2(V+A) | p(V+A) | Chi2(V+A+P)  | p(V+A+P) '; 
out = strvcat(out, sprintf('--------------------------------------------------------------------------------------------------------'));    
clear OutputValue;OutputValue=[chi2_Vel chi2P_Vel chi2_VelAcc chi2P_VelAcc chi2_VelAccPos chi2P_VelAccPos];
% clear OutputValue;OutputValue=[RChi2_Vel chi2P_Vel RChi2_VelAcc chi2P_VelAcc RChi2_VelAccPos chi2P_VelAccPos];
out=strvcat(out, sprintf('                   %6.3f  | %6.3f |  %6.3f   | %6.3f |  %6.3f      | %6.3f    ', OutputValue));
figure(FigureIndex);axes('position',[0.26 0.59 0.25 0.25]); set(gca,'box','off','visible','off');
text(-0.08,1,out,'fontsize',8,'fontname','courier','horizontalalignment','left','verticalalignment','top');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do AIC test
AIC_1vs2=size()*log(err_VelAcc/err_Vel)+2*(length(vect_VelAcc)-length(vect_Vel));
AIC_1vs2=234*log(err_VelAcc/err_Vel)+2*(length(vect_VelAcc)-length(vect_Vel));

234*log(374.16/536.83)+2*(length(vect_VelAcc)-length(vect_Vel))
234*log(373.19/374.16)+2*(length(vect_VelAccPos)-length(vect_VelAcc))

AIC = - 2*log L + k * edf,

%    I've been doing some self-study of forecasting and have a question
%     about calculating the AIC. I have a time series of data and wish to
%     determine the linear model of order p which is most appropriate for
%     the data. Suppose the time series is {r_t} with t=1,2,...,T. I am
%     fitting a linear model of the form: 
%     r_t = x_0 + x_1 r_(t-1) + ... + x_p r_(t-p) + e_t 
%     I have seen several different definitions of the AIC, most commonly
%     AIC(k) = Log(sigmahat_k^2) + 2k/T
%     I think my question is on calculating sigmahat_k^2. Do I have to use
%     least squares to estimate the parameters x_0, x_1, ..., x_p, then
%     calculate e_t for t=k+1, k+2, ..., T, and then find the sum of squares
%     of the e_t, and then divide by T-p? Or is there an easier way to do
%     this?

% You under stand the problem correctly. Sigmahat_p^2 is the mean square
% for errors: residual sum of squares divided by the number of
% observations (N) minus the number of regression paramters (p+1)). It's
% the average of the squared residuals, e_t, but using N-p-1 as the
% divisor instead of N. Then plot AIC(p) against p. It should level off
% and pick the value of p where it starts to level off. Additional
% parmeters do not imporve the fit.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% order in which the directions are plotted
plot_col = [1 1 1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 6 6 6 7 7 7 8 8 8];
plot_row = [5 4 3 2 1 4 3 2 4 3 2 4 3 2 4 3 2 4 3 2 4 3 2 4 3 2];

x_start = [0, 0];
x_stop =  [2, 2];
y_marker=[0,1.1*max(max(spacetime_data))];

xscale = [0.8 0.77 0.65 0.52 0.5 0.52 0.65 0.77];   
yscale = [0.35 0.5 0.65 0.5 0.35 0.2 0.05 0.2];

%x_time=[1:size(PSTH_tempPlane,2)]*0.02-0.5;
clear x_time;x_time=[1:size(spacetime_data,2)]*timestep;

UniAzi=[0:45:315];
for i=1:8
    i
    axes('position',[xscale(i) yscale(i) 0.16 0.10]);
%     bar(x_time,PSTH_tempPlane(i,:));hold on;%bar(x_time, count_y{i,j,k});    hold on;    
    plot(x_time,spacefit_Vel(i,:),'r','LineWidth',2);hold on;    
    plot(x_time,spacefit_VelAcc(i,:), 'g', 'LineWidth', 2);hold on;
    plot(x_time,spacefit_VelAccPos(i,:),'c','LineWidth',2);hold on;
    
     xlim([0,2.5]);    ylim([0,1.2*max(max(spacetime_data))]);
    
    set(gca, 'xtick', [] );set(gca, 'XTickMode','manual'); 
    set(gca, 'xtick',[0:0.5:2.5]);
    set( gca, 'xticklabel', '0|0.5|1|1.5|2|2.5');
    set( gca, 'yticklabel', ' ' );
    title({['azi = ' num2str(UniAzi(i))]}, 'FontSize', 8);
end  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save the figure 
OutputPath=['C:\Aihua\z_TempOutputs\'];
OutputPath=['Z:\Users\Aihua\CurveFittingAnalysis\outputs\'];
figure(FigureIndex); 
set(gcf, 'PaperOrientation', 'portrait');
saveas(gcf,[OutputPath 'FakeData_Pos.png'],'png');
close(FigureIndex);






