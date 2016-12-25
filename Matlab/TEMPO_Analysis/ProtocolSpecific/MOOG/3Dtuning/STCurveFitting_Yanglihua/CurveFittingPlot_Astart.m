function CurveFit=CurveFittingPlot_Astart(FILE,Plane,Azi_3D,Ele_3D,Azi_temp,Ele_temp,Step,x_time,x_stop,StepMatrix,StdMatrix,SpikeCount_Trial,p_peak,Value_peak,TimeIndex_peak,p_trough,Value_trough,TimeIndex_trough,Sti)

% Step=100;
Azi_3D=[0:45:315];Ele_3D=[-90:45:90];
[AziGrid,EleGrid]=meshgrid(Azi_3D,Ele_3D); 
% x_timeL=[0:-Step*0.001:-0.5];x_timeL=fliplr(x_timeL);x_timeL=x_timeL(1:end-1);
% % x_timeR=[0:Step*0.001:2.1];%
% x_timeR=[0:Step*0.001:2.5];
% x_time=[x_timeL x_timeR];
t=[0:Step*0.001:x_time(end)];
clear StartIndex; StartIndex=find(x_time==0);
clear XAzi Ytime;[XAzi,Ytime] = meshgrid([0:45:360],[0:Step*0.001:max(t)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get the space time data for different plane (only for vestibular condition)
for i=1:length(Azi_temp) 
    clear RowIndex ColIndex; [RowIndex, ColIndex]=find(AziGrid==Azi_temp(i) & EleGrid==Ele_temp(i));
    PSTH_tempPlane(i,:)=StepMatrix{Sti}(RowIndex,ColIndex,:);  %PSTH_tempPlane(i,:)=StepMatrix{1}(RowIndex,ColIndex,:);  
    p_peak_tempPlane(i,1)=p_peak{Sti}(RowIndex,ColIndex);
    p_trough_tempPlane(i,1)=p_trough{Sti}(RowIndex,ColIndex);
    PeakValue_tempPlane(i,1)=Value_peak{Sti}(RowIndex,ColIndex);
    TroughValue_tempPlane(i,1)=Value_trough{Sti}(RowIndex,ColIndex);
    PeakTimeIndex(i,1)=TimeIndex_peak{Sti}(RowIndex,ColIndex);
    if x_time(PeakTimeIndex)>2
        p_peak_tempPlane(i,1)=NaN;
    else
    end
    TroughTimeIndex(i,1)=TimeIndex_trough{Sti}(RowIndex,ColIndex); 
    if x_time(TroughTimeIndex(i,1))>2
        p_trough_tempPlane(i,1)=NaN;
    else
    end    
    clear tempdata; tempdata=[SpikeCount_Trial{Sti,RowIndex,ColIndex}]';
 
    if i>1
        tempdata0=zeros(size(spacetime_data_trial,1),size(spacetime_data_trial,3));
        tempdata0(1:size(tempdata,1),:)=tempdata(:,StartIndex:end);
        spacetime_data_trial(:,i,:)=tempdata0;
    else
        spacetime_data_trial(:,i,:)=tempdata(:,StartIndex:end);
    end    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Check the temporal modulation 
clear Index_Sig; Index_Sig=find(p_peak_tempPlane<0.01);
k=0;
for i=1:length(Index_Sig)-1
    for j=i+1:length(Index_Sig)
        k=k+1;
        DiffDir_Peak(k,1)=Angle3D_paired(Azi_temp(Index_Sig(i)),Azi_temp(Index_Sig(j)), Ele_temp(Index_Sig(i)),Ele_temp(Index_Sig(j)));
    end
end
if length(Index_Sig)>=2 & min(DiffDir_Peak)<=50
    Modulation_Peak=1;
else
    Modulation_Peak=0;
end

clear Index_Sig; Index_Sig=find(p_trough_tempPlane<0.01);
k=0;
for i=1:length(Index_Sig)-1
    for j=i+1:length(Index_Sig)
        k=k+1;
        DiffDir_Trough(k,1)=Angle3D_paired(Azi_temp(Index_Sig(i)),Azi_temp(Index_Sig(j)), Ele_temp(Index_Sig(i)),Ele_temp(Index_Sig(j)));
    end
end
if length(Index_Sig)>=2 & min(DiffDir_Trough)<=50
    Modulation_Trough=1;
else
    Modulation_Trough=0;
end
% if Modulation_Peak==0 & Modulation_Trough==0
%     CurveFit=[];
%     return;
% else
% end  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot the original space-time figure; 
for i=1:size(spacetime_data_trial,1)
    clear tempdata; tempdata=squeeze(spacetime_data_trial(i,:,:));
    tempdata(9,:)=tempdata(1,:);
    tempdata=reshape(tempdata,1,size(tempdata,1)*size(tempdata,2));
    spacetime_data_re(i,:)=tempdata;%    spacetime_data_re(i,:)=tempdata*1000/Step;
    Direction_re(i,:)=reshape(XAzi,1,size(tempdata,1)*size(tempdata,2));
    time_re(i,:)=reshape(Ytime,1,size(tempdata,1)*size(tempdata,2));
end
spacetime_data=reshape(mean(spacetime_data_re),length(Azi_temp)+1,length(t));
%%%%%%%%%%%%%%%%%
% Do the Two-way ANOVA to see whether the space-time data has a significant structure
clear currentdata; currentdata=[];
for i=1:size(spacetime_data_trial,2)
    clear tempdata;tempdata=squeeze(spacetime_data_trial(:,i,:));
    currentdata=cat(1,currentdata,tempdata);
end
[p_anova_2way,tbl,stats] = anova2(currentdata,size(spacetime_data_trial,1),'off');

if p_anova_2way(1)>0.001 | p_anova_2way(2)>0.001 | p_anova_2way(3)>0.001
    CurveFit=[];
    return;
end  

FigureIndex=2; 
figure(FigureIndex);set(FigureIndex,'Position', [50,50 1200,800], 'Name', 'CurveFitting');orient landscape; 
text(-0.1,1.06,[FILE '  ' Plane]); axis off;    
clear XAzi Ytime;[XAzi,Ytime] = meshgrid([0:45:360],[0:Step*0.001:max(t)]);
figure(FigureIndex);axes('position',[0.05 0.77 0.17 0.15]);%subplot('position', [0.1 0.7 0.22 0.22]);
contourf(XAzi,Ytime,spacetime_data');
% caxis([0 120]);
colorbar;    
set(gca, 'xtick', [] );set(gca, 'XTickMode','manual'); 
set(gca, 'xtick',[0:90:360]);set(gca, 'xticklabel','0|90|180|270|360');
set(gca, 'ytick', [] ); set(gca, 'YTickMode','manual');
set(gca, 'ytick',[0:0.5:max(t)]);
set(gca, 'yticklabel','0|0.5|1|1.5|2|2.5');% set(gca, 'ydir','reverse');
xlabel('Azimuth, X (deg)');  ylabel('Time (sec)');title('raw data');

out=[];
out = 'Two-way ANOVA: p(Time) | p(Space) | p(Space*Time) '; 
out = strvcat(out, sprintf('--------------------------------------------------------------------------------------------------------'));    
clear OutputValue;OutputValue=p_anova_2way;
out=strvcat(out, sprintf('              %7.3f  | %7.3f |  %7.3f  ', OutputValue));
figure(FigureIndex);axes('position',[0.26 0.74 0.25 0.25]); set(gca,'box','off','visible','off');
text(-0.08,1,out,'fontsize',8,'fontname','courier','horizontalalignment','left','verticalalignment','top');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Do the curve fitting
%Model1: Acceleration only
allow_negative=0;
global model_use
model_use=0;
[spacefit_Acc,vect_Acc,r_squared_Acc,CI_Acc,cor_Acc] = MSF_Acc_fit(spacetime_data,Step*0.001,allow_negative);

figure(FigureIndex);axes('position',[0.05 0.53 0.17 0.15]);
contourf(XAzi,Ytime,spacefit_Acc',10);%caxis([0 120]);
colorbar;    
set(gca, 'xtick', [] );set(gca, 'XTickMode','manual'); 
set(gca, 'xtick',[0:90:360]);set(gca, 'xticklabel','0|90|180|270|360');
set(gca, 'ytick', [] ); set(gca, 'YTickMode','manual');
set(gca, 'ytick',[0:0.5:max(t)]);
set(gca, 'yticklabel','0|0.5|1|1.5|2|2.5');% set(gca, 'ydir','reverse');
ylabel('Time (sec)');title('model: Acc  ');%xlabel('Azimuth, X (deg)');  

error_surf = spacetime_data - spacefit_Acc;
err_Acc = cosnlin_err(vect_Acc);
figure(FigureIndex);axes('position',[0.26 0.53 0.17 0.15]);%subplot('position', [0.1 0.7 0.22 0.22]);
contourf(XAzi,Ytime,error_surf'); 
% caxis([0 120]);
colorbar; 
set(gca, 'xtick', [] );set(gca, 'XTickMode','manual'); 
set(gca, 'xtick',[0:90:360]);set(gca, 'xticklabel','0|90|180|270|360');
set(gca, 'ytick', [] ); set(gca, 'YTickMode','manual');
set(gca, 'ytick',[0:0.5:max(t)]);
set(gca, 'yticklabel','0|0.5|1|1.5|2|2.5');% set(gca, 'ydir','reverse');
%axis off;
title({[ 'Err: ' num2str(err_Acc, '%0.2f') ]},  'FontSize', 10);

% Do chi-square goodness of fit test
%model 1
global xdata tdata 
xtdata = [xdata;tdata];
model_use=1;
[chi2_Acc, chi2P_Acc] = Chi2_Test_3D(xtdata, spacetime_data_re, 'funccosnlin', vect_Acc, length(vect_Acc));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Model 2: Velocity + Acceleration
model_use=2;
[spacefit_VelAcc,vect_VelAcc,r_squared_VelAcc,CI_VelAcc,cor_VelAcc] = MSF_VelAcc_fit(spacetime_data, vect_Acc,Step*0.001,allow_negative);

figure(FigureIndex);axes('position',[0.05 0.31 0.17 0.15]);
contourf(XAzi,Ytime,spacefit_VelAcc');%caxis([0 120]);
colorbar;    
set(gca, 'xtick', [] );set(gca, 'XTickMode','manual'); 
set(gca, 'xtick',[0:90:360]);set(gca, 'xticklabel','0|90|180|270|360');
set(gca, 'ytick', [] ); set(gca, 'YTickMode','manual');
set(gca, 'ytick',[0:0.5:max(t)]);
set(gca, 'yticklabel','0|0.5|1|1.5|2|2.5');% set(gca, 'ydir','reverse');
ylabel('Time (sec)');title('model: Vel + Acc ');%xlabel('Azimuth, X (deg)');  

error_surf = spacetime_data - spacefit_VelAcc;
err_VelAcc = cosnlin_err(vect_VelAcc);
figure(FigureIndex);axes('position',[0.26 0.31 0.17 0.15]);
contourf(XAzi,Ytime,error_surf'); %caxis([0 120]);
colorbar;
set(gca, 'xtick', [] );set(gca, 'XTickMode','manual'); 
set(gca, 'xtick',[0:90:360]);set(gca, 'xticklabel','0|90|180|270|360');
set(gca, 'ytick', [] ); set(gca, 'YTickMode','manual');
set(gca, 'ytick',[0:0.5:max(t)]);
set(gca, 'yticklabel','0|0.5|1|1.5|2|2.5');% set(gca, 'ydir','reverse');
% axis off;
title({[ 'Err: ' num2str(err_VelAcc, '%0.2f') ]},  'FontSize', 10);

% Do chi-square goodness of fit test 
[chi2_VelAcc, chi2P_VelAcc] = Chi2_Test_3D(xtdata, spacetime_data_re, 'funccosnlin', vect_VelAcc, length(vect_VelAcc));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Model 3: Velocity + Acceleration + Position
model_use=3;
[spacefit_VelAccPos,vect_VelAccPos,r_squared_VelAccPos,CI_VelAccPos,cor_VelAccPos] = MSF_VelAccPos_fit(spacetime_data, vect_VelAcc,Step*0.001,allow_negative);

figure(FigureIndex);axes('position',[0.05 0.09 0.17 0.15]);%subplot('position', [0.1 0.7 0.22 0.22]);
contourf(XAzi,Ytime,spacefit_VelAccPos',10);%caxis([0 120]);
colorbar;    
set(gca, 'xtick', [] );set(gca, 'XTickMode','manual'); 
set(gca, 'xtick',[0:90:360]);set(gca, 'xticklabel','0|90|180|270|360');
set(gca, 'ytick', [] ); set(gca, 'YTickMode','manual');
set(gca, 'ytick',[0:0.5:max(t)]);
set(gca, 'yticklabel','0|0.5|1|1.5|2|2.5');% set(gca, 'ydir','reverse');
xlabel('Azimuth, X (deg)'); 
ylabel('Time (sec)');title('model: Vel + Acc + Pos');

error_surf = spacetime_data - spacefit_VelAccPos;
err_VelAccPos = cosnlin_err(vect_VelAccPos);
figure(FigureIndex);axes('position',[0.26 0.09 0.17 0.15]);
contourf(XAzi,Ytime,error_surf'); %caxis([0 120]);
colorbar;
set(gca, 'xtick', [] );set(gca, 'XTickMode','manual'); 
set(gca, 'xtick',[0:90:360]);set(gca, 'xticklabel','0|90|180|270|360');
set(gca, 'ytick', [] ); set(gca, 'YTickMode','manual');
set(gca, 'ytick',[0:0.5:max(t)]);
set(gca, 'yticklabel','0|0.5|1|1.5|2|2.5');% set(gca, 'ydir','reverse');
% axis off; 
title({[ 'Err: ' num2str(err_VelAccPos, '%0.2f') ]},  'FontSize', 10);% axis off;

% mode 3
[chi2_VelAccPos, chi2P_VelAccPos] = Chi2_Test_3D(xtdata, spacetime_data_re, 'funccosnlin', vect_VelAccPos, length(vect_VelAccPos));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out=[];
out = 'CurveFitting:    R0 |  Amp  |  n  | muAzi |  muT  | sigmaT |  DC2  | wVel  |ThetaAcc|  wAcc | ThetaPos | wPos'; 
out = strvcat(out, sprintf('------------------------------------------------------------------------------------------------------------'));    
clear OutputValue;OutputValue=vect_Acc;
OutputValue(4)=OutputValue(4)*180/pi;
out=strvcat(out, sprintf('Acc Model:     %4.1f | %5.1f |%4.1f |%6.1f | %4.3f | %6.3f |%6.3f', OutputValue));

clear OutputValue;OutputValue=vect_VelAcc;
OutputValue(4)=OutputValue(4)*180/pi;
OutputValue(9)=OutputValue(9)*180/pi;
out=strvcat(out, sprintf('Vel+Acc Model: %4.1f | %5.1f |%4.1f |%6.1f | %4.3f | %6.3f |%6.3f |%6.3f | %6.1f  ', OutputValue));

clear OutputValue; OutputValue=[vect_VelAccPos(1:9) 1-vect_VelAccPos(8) vect_VelAccPos(10:11)];
OutputValue(4)=OutputValue(4)*180/pi;
OutputValue(9)=OutputValue(9)*180/pi;
OutputValue(11)=OutputValue(11)*180/pi;
out=strvcat(out, sprintf('Vel+Acc+Pos  : %4.1f | %5.1f |%4.1f |%6.1f | %4.3f | %6.3f |%6.3f |%6.3f | %6.1f |%6.3f | %6.1f   |%6.3f', OutputValue)); 
figure(FigureIndex);axes('position',[0.26 0.69 0.25 0.25]); set(gca,'box','off','visible','off');
text(-0.08,1,out,'fontsize',8,'fontname','courier','horizontalalignment','left','verticalalignment','top');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Performing regression to get R^2 and p-values for F-test
Data_Raw=reshape(spacetime_data,size(spacetime_data,1)*size(spacetime_data,2),1);
Data_Acc=reshape(spacefit_Acc,size(spacefit_Acc,1)*size(spacefit_Acc,2),1);
Data_VelAcc=reshape(spacefit_VelAcc,size(spacefit_VelAcc,1)*size(spacefit_VelAcc,2),1);
Data_VelAccPos=reshape(spacefit_VelAccPos,size(spacefit_VelAccPos,1)*size(spacefit_VelAccPos,2),1);

clear X1;X1 = [ones(size(Data_Acc,1),1) Data_Acc];%X1 = [ones(808,1) Data_Vel];% y_fit = [ones(length(y_fit),1) y_fit];
[b_Acc,bint,r_Acc,rint,stats_Acc] = regress(Data_Raw,X1,0.05);
clear X1;X1 = [ones(size(Data_VelAcc,1),1) Data_VelAcc];
[b_VelAcc,bint,r_VelAcc,rint,stats_VelAcc] = regress(Data_Raw,X1,0.05);
clear X1;X1 = [ones(size(Data_VelAccPos,1),1) Data_VelAccPos];
[b_VelAccPos,bint,r_VelAccPos,rint,stats_VelAccPos] = regress(Data_Raw,X1,0.05);

Ftest_1vs2=[(err_Acc-err_VelAcc)/(length(vect_VelAcc)-length(vect_Acc))]/[err_VelAcc/(size(spacetime_data,1)*size(spacetime_data,2)-length(vect_VelAcc))];
p_1vs2=1-fcdf(Ftest_1vs2,length(vect_VelAcc)-length(vect_Acc),size(spacetime_data,1)*size(spacetime_data,2)-length(vect_VelAcc));

Ftest_2vs3=[(err_VelAcc-err_VelAccPos)/(length(vect_VelAccPos)-length(vect_VelAcc))]/[err_VelAccPos/(size(spacetime_data,1)*size(spacetime_data,2)-length(vect_VelAccPos))];
p_2vs3=1-fcdf(Ftest_2vs3,length(vect_VelAccPos)-length(vect_VelAcc),size(spacetime_data,1)*size(spacetime_data,2)-length(vect_VelAccPos));

% do AIC test
AIC_1vs2=size(spacetime_data,1)*size(spacetime_data,2)*log(err_VelAcc/err_Acc)+2*(length(vect_VelAcc)-length(vect_Acc));
AIC_2vs3=size(spacetime_data,1)*size(spacetime_data,2)*log(err_VelAccPos/err_VelAcc)+2*(length(vect_VelAccPos)-length(vect_VelAcc));
AIC_1vs3=size(spacetime_data,1)*size(spacetime_data,2)*log(err_VelAccPos/err_Acc)+2*(length(vect_VelAccPos)-length(vect_Acc));

out=[];
% out = 'Goodness of fitting:  r2 (Vel) | r2 (VelAcc)| r2 (VelAccPos) | Ftest(1vs2) | p(1vs2)  | Ftest(2vs3) | p(2vs3)'; 
out = 'Goodness of fit: r2(Acc) | r2(V+A)| r2(V+A+P) | F(1vs2) | p(1vs2)  | F(2vs3) | p(2vs3) | AIC(1vs2) | AIC(2vs3)'; 
out = strvcat(out, sprintf('--------------------------------------------------------------------------------------------------------'));    
clear OutputValue;OutputValue=[r_squared_Acc r_squared_VelAcc r_squared_VelAccPos Ftest_1vs2 p_1vs2 Ftest_2vs3 p_2vs3 AIC_1vs2 AIC_2vs3];
out=strvcat(out, sprintf('                  %6.3f | %6.3f |  %6.3f   | %6.3f  |  %6.3f  |  %6.3f |  %6.3f | %9.3f |  %6.3f ', OutputValue));
figure(FigureIndex);axes('position',[0.26 0.60 0.25 0.25]); set(gca,'box','off','visible','off');
text(-0.08,1,out,'fontsize',8,'fontname','courier','horizontalalignment','left','verticalalignment','top');

out=[];
out = 'Goodness of fit: Chi2(A) |  p(A) | Chi2(VA)|  p(VA) | Chi2(VAP) | p(VAP) '; 
out = strvcat(out, sprintf('---------------------------------------------------------------------------------'));    
clear OutputValue;OutputValue=[chi2_Acc chi2P_Acc chi2_VelAcc chi2P_VelAcc chi2_VelAccPos chi2P_VelAccPos];
out=strvcat(out, sprintf('             %9.1f   | %5.3f | %7.1f |  %5.3f | %8.1f  | %6.3f   ', OutputValue));
figure(FigureIndex);axes('position',[0.26 0.55 0.25 0.25]); set(gca,'box','off','visible','off');
text(-0.08,1,out,'fontsize',8,'fontname','courier','horizontalalignment','left','verticalalignment','top');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% order in which the directions are plotted
plot_col = [1 1 1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 6 6 6 7 7 7 8 8 8];
plot_row = [5 4 3 2 1 4 3 2 4 3 2 4 3 2 4 3 2 4 3 2 4 3 2 4 3 2];

x_start = [0, 0];
% x_stop =  [2, 2];
y_marker=[0,1.1*max(max(spacetime_data))];

xscale = [0.8 0.77 0.68 0.52 0.5 0.52 0.68 0.77];  
yscale = [0.34 0.49 0.64 0.49 0.34 0.19 0.04 0.19];

% xscale = [0.8 0.77 0.65 0.52 0.5 0.52 0.65 0.77];  
% yscale = [0.35 0.5 0.65 0.5 0.35 0.2 0.05 0.2];

%get the curve fitting results for each direction
clear model_Acc; model_Acc=zeros(length(Azi_temp),length(x_time));
model_Acc(:,1:StartIndex-1) = NaN; model_Acc(:,StartIndex:length(x_time))=spacefit_Acc(1:8,:);
clear model_VelAcc; model_VelAcc=zeros(length(Azi_temp),length(x_time));
model_VelAcc(:,1:StartIndex-1) = NaN; model_VelAcc(:,StartIndex:length(x_time))=spacefit_VelAcc(1:8,:);
clear model_VelAccPos; model_VelAccPos=zeros(length(Azi_temp),length(x_time));
model_VelAccPos(:,1:StartIndex-1) = NaN; model_VelAccPos(:,StartIndex:length(x_time))=spacefit_VelAccPos(1:8,:);

%Plot the PSTH, superimposed with curve fitting
for i=1:length(Azi_temp)
    i;
    axes('position',[xscale(i) yscale(i) 0.16 0.10]);
    bar(x_time,PSTH_tempPlane(i,:),1.0);hold on;%bar(x_time, count_y{i,j,k});    hold on;
    
     clear temp_p_peak; temp_p_peak=p_peak_tempPlane(i,1);
     if temp_p_peak <0.01 & x_time(PeakTimeIndex(i,1))<2            
         plot(x_time(PeakTimeIndex(i,1)), PSTH_tempPlane(i,PeakTimeIndex(i,1)),'ro','LineWidth',2.0);hold on;
     end

     clear temp_p_trough; temp_p_trough=p_trough_tempPlane(i,1);
     if temp_p_trough <0.01 & x_time(TroughTimeIndex(i,1))<2           
         plot(x_time(TroughTimeIndex(i,1)), PSTH_tempPlane(i,TroughTimeIndex(i,1)),'go','LineWidth',2.0);hold on;
     end
   
    plot(x_time,model_Acc(i,:),'r','LineWidth',2); hold on;   
    plot(x_time,model_VelAcc(i,:), 'g', 'LineWidth', 2);hold on;
    plot(x_time,model_VelAccPos(i,:),'c','LineWidth',2);hold on;
    
    plot(x_start, y_marker, 'k-','LineWidth',2);hold on;
    plot(x_stop,  y_marker, 'k-','LineWidth',2);hold on;
    xlim([-0.5,2.5]);    ylim([0,1.2*max(max(spacetime_data))]);
    
    set(gca, 'xtick', [] );set(gca, 'XTickMode','manual'); 
    set(gca, 'xtick',[-0.5:0.5:2.5]);
    set(gca, 'xticklabel', '-0.5|0|0.5|1|1.5|2|2.5');
    set(gca, 'yticklabel', ' ' );
    title({['azi=' num2str(Azi_temp(i)) '; Ele=' num2str(Ele_temp(i))]}, 'FontSize', 8);
end  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output some data for further analysis
CurveFit.spacetime_data=spacetime_data;
CurveFit.Modulation_Peak=Modulation_Peak;
CurveFit.Modulation_Trough=Modulation_Trough;

CurveFit.PeakValue=PeakValue_tempPlane;
CurveFit.TroughValue=TroughValue_tempPlane;
CurveFit.p_peak=p_peak_tempPlane;
CurveFit.p_trough=p_trough_tempPlane;

CurveFit.p_anova_2way=p_anova_2way;

CurveFit.spacefit_Acc=spacefit_Acc;
CurveFit.spacefit_VelAcc=spacefit_VelAcc;
CurveFit.spacetime_VelAccPos=spacefit_VelAccPos;
CurveFit.vect_Acc=vect_Acc;
CurveFit.vect_VelAcc=vect_VelAcc;
CurveFit.vect_VelAccPos=vect_VelAccPos;
CurveFit.CI_Acc=CI_Acc;
CurveFit.CI_VelAcc=CI_VelAcc;
CurveFit.CI_VelAccPos=CI_VelAccPos;

CurveFit.cor_Acc=cor_Acc;
CurveFit.cor_VelAcc=cor_VelAcc;
CurveFit.cor_VelAccPos=cor_VelAccPos;

CurveFit.stats_Acc=stats_Acc;
CurveFit.stats_VelAcc=stats_VelAcc;
CurveFit.stats_VelAccPos=stats_VelAccPos;
CurveFit.err_Acc=err_Acc;
CurveFit.err_VelAcc=err_VelAcc;
CurveFit.err_VelAccPos=err_VelAccPos;
CurveFit.r_squared_Acc=r_squared_Acc;
CurveFit.r_squared_VelAcc=r_squared_VelAcc;
CurveFit.r_squared_VelAccPos=r_squared_VelAccPos;
CurveFit.Ftest_1vs2=Ftest_1vs2;
CurveFit.p_1vs2=p_1vs2;
CurveFit.Ftest_2vs3=Ftest_2vs3;
CurveFit.p_2vs3=p_2vs3;

CurveFit.AIC_1vs2=AIC_1vs2;
CurveFit.AIC_2vs3=AIC_2vs3;

CurveFit.chi2_Acc=chi2_Acc;% CurveFit.Rchi2_Vel=RChi2_Vel;
CurveFit.chi2P_Acc=chi2P_Acc;% CurveFit.dof_Vel=dof_Vel;
CurveFit.chi2_VelAcc=chi2_VelAcc;% CurveFit.Rchi2_VelAcc=RChi2_VelAcc;
CurveFit.chi2P_VelAcc=chi2P_VelAcc;%CurveFit.dof_VelAcc=dof_VelAcc;
CurveFit.chi2_VelAccPos=chi2_VelAccPos;%CurveFit.Rchi2_VelAccPos=RChi2_VelAccPos;
CurveFit.chi2P_VelAccPos=chi2P_VelAccPos;%CurveFit.dof_VelAccPos=dof_VelAccPos;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Save the figure 
% OutputPath=['C:\Aihua\z_TempOutputs\'];
% figure(FigureIndex); 
% set(gcf, 'PaperOrientation', 'portrait');
% saveas(gcf,[OutputPath FILE(1:end-4) '_CurveFitting.png'],'png');
% close(FigureIndex);
% 
% %Save the Data 
% SaveFileName=[OutputPath FILE(1:end-4) '_CurveFit'];
% save(SaveFileName,'CurveFit'); clear SaveFileName;

