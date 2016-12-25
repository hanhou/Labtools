function CurveFit=CurveFitting_Plot_Vel_NDC2(FILE,Plane,Azi_3D,Ele_3D,Azi_temp,Ele_temp,Step,StepMatrix,StdMatrix,SpikeCount_Trial,p_peak,Value_peak,TimeIndex_peak,p_trough,Value_trough,TimeIndex_trough)

% Step=100;
Azi_3D=[0:45:315];Ele_3D=[-90:45:90];
[AziGrid,EleGrid]=meshgrid(Azi_3D,Ele_3D); 
x_timeL=[0:-Step*0.001:-0.5];x_timeL=fliplr(x_timeL);x_timeL=x_timeL(1:end-1);
% x_timeR=[0:Step*0.001:2.1];%
x_timeR=[0:Step*0.001:2.5];
x_time=[x_timeL x_timeR];
t=[0:Step*0.001:x_time(end)];
clear StartIndex; StartIndex=find(x_time==0);
clear XAzi Ytime;[XAzi,Ytime] = meshgrid([0:45:360],[0:Step*0.001:max(t)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get the space time data for different plane (only for vestibular condition)
for i=1:length(Azi_temp) 
    clear RowIndex ColIndex; [RowIndex, ColIndex]=find(AziGrid==Azi_temp(i) & EleGrid==Ele_temp(i));
    PSTH_tempPlane(i,:)=StepMatrix{1}(RowIndex,ColIndex,:);%PSTH_tempPlane(i,:)=StepMatrix{1}(RowIndex,ColIndex,:)*1000/Step;  
    p_peak_tempPlane(i,1)=p_peak{1}(RowIndex,ColIndex);
    p_trough_tempPlane(i,1)=p_trough{1}(RowIndex,ColIndex);
    PeakValue_tempPlane(i,1)=Value_peak{1}(RowIndex,ColIndex);
    TroughValue_tempPlane(i,1)=Value_trough{1}(RowIndex,ColIndex);
    PeakTimeIndex(i,1)=TimeIndex_peak{1}(RowIndex,ColIndex);
    if x_time(PeakTimeIndex)>2
        p_peak_tempPlane(i,1)=NaN;
    else
    end
    TroughTimeIndex(i,1)=TimeIndex_trough{1}(RowIndex,ColIndex); 
    if x_time(TroughTimeIndex(i,1))>2
        p_trough_tempPlane(i,1)=NaN;
    else
    end    
    clear tempdata; tempdata=[SpikeCount_Trial{1,RowIndex,ColIndex}]';%     SpikeCount_Trial_tempPlane{i}=[SpikeCount_Trial{1,RowIndex,ColIndex}]';
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
if Modulation_Peak==0 & Modulation_Trough==0
    CurveFit=[];
    return;
else
end  

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
figure(FigureIndex);set(FigureIndex,'Position', [50,100 1200,800], 'Name', 'CurveFitting');orient landscape; 
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
out=strvcat(out, sprintf('               %7.3f  | %7.3f |  %7.3f  ', OutputValue));
figure(FigureIndex);axes('position',[0.26 0.59 0.25 0.25]); set(gca,'box','off','visible','off');
text(-0.08,1,out,'fontsize',8,'fontname','courier','horizontalalignment','left','verticalalignment','top');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Do the curve fitting 
% Vel Model (without adding DC2)
allow_negative=0;
global model_use
model_use=1;
[spacefit_Vel_NDC2,vect_Vel_NDC2,r_squared_Vel_NDC2,CI_Vel_NDC2,cor_Vel_NDC2] = MSF_Vel_fit_NDC2(spacetime_data,Step*0.001,allow_negative);

figure(FigureIndex);axes('position',[0.05 0.53 0.17 0.15]);
contourf(XAzi,Ytime,spacefit_Vel_NDC2',10);%caxis([0 120]);
colorbar;    
set(gca, 'xtick', [] );set(gca, 'XTickMode','manual'); 
set(gca, 'xtick',[0:90:360]);set(gca, 'xticklabel','0|90|180|270|360');
set(gca, 'ytick', [] ); set(gca, 'YTickMode','manual');
set(gca, 'ytick',[0:0.5:max(t)]);
set(gca, 'yticklabel','0|0.5|1|1.5|2|2.5');% set(gca, 'ydir','reverse');
ylabel('Time (sec)');title('model: Vel (Without DC2)');%xlabel('Azimuth, X (deg)');  

error_surf = spacetime_data - spacefit_Vel_NDC2;
err_Vel_NDC2 = cosnlin_err_NDC2(vect_Vel_NDC2);
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
title({[ 'Err: ' num2str(err_Vel_NDC2, '%0.2f') ]},  'FontSize', 10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Vel Model (adding DC2)
[spacefit_Vel,vect_Vel,r_squared_Vel,CI_Vel,cor_Vel] = MSF_Vel_fit(spacetime_data,Step*0.001,allow_negative);

figure(FigureIndex);axes('position',[0.05 0.31 0.17 0.15]);
contourf(XAzi,Ytime,spacefit_Vel');%caxis([0 120]);
colorbar;    
set(gca, 'xtick', [] );set(gca, 'XTickMode','manual'); 
set(gca, 'xtick',[0:90:360]);set(gca, 'xticklabel','0|90|180|270|360');
set(gca, 'ytick', [] ); set(gca, 'YTickMode','manual');
set(gca, 'ytick',[0:0.5:max(t)]);
set(gca, 'yticklabel','0|0.5|1|1.5|2|2.5');% set(gca, 'ydir','reverse');
ylabel('Time (sec)');title('model: Vel (with DC2)');%xlabel('Azimuth, X (deg)');  

error_surf = spacetime_data - spacefit_Vel;
err_Vel = cosnlin_err(vect_Vel);
figure(FigureIndex);axes('position',[0.26 0.31 0.17 0.15]);
contourf(XAzi,Ytime,error_surf'); %caxis([0 120]);
colorbar;
set(gca, 'xtick', [] );set(gca, 'XTickMode','manual'); 
set(gca, 'xtick',[0:90:360]);set(gca, 'xticklabel','0|90|180|270|360');
set(gca, 'ytick', [] ); set(gca, 'YTickMode','manual');
set(gca, 'ytick',[0:0.5:max(t)]);
set(gca, 'yticklabel','0|0.5|1|1.5|2|2.5');% set(gca, 'ydir','reverse');
% axis off;
title({[ 'Err: ' num2str(err_Vel, '%0.2f') ]},  'FontSize', 10);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out=[];
out = 'CurveFitting:          R0 |  Amp  |  n  | muAzi |  muT  | sigmaT |  DC2  '; 
out = strvcat(out, sprintf('------------------------------------------------------------------------------------------------------------'));    
clear OutputValue;OutputValue=vect_Vel_NDC2;
OutputValue(4)=OutputValue(4)*180/pi;
out=strvcat(out, sprintf('Vel Model(No DC2):   %4.1f | %5.1f |%4.1f |%6.1f | %4.3f | %6.3f ', OutputValue));

clear OutputValue;OutputValue=vect_Vel;
OutputValue(4)=OutputValue(4)*180/pi;
out=strvcat(out, sprintf('Vel Model(With DC2): %4.1f | %5.1f |%4.1f |%6.1f | %4.3f | %6.3f |%6.3f ', OutputValue));

figure(FigureIndex);axes('position',[0.26 0.74 0.25 0.25]); set(gca,'box','off','visible','off');
text(-0.08,1,out,'fontsize',8,'fontname','courier','horizontalalignment','left','verticalalignment','top');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Performing regression to get R^2 and p-values for F-test
Data_Raw=reshape(spacetime_data,size(spacetime_data,1)*size(spacetime_data,2),1);
Data_Vel=reshape(spacefit_Vel,size(spacefit_Vel,1)*size(spacefit_Vel,2),1);
Data_Vel_NDC2=reshape(spacefit_Vel_NDC2,size(spacefit_Vel,1)*size(spacefit_Vel,2),1);


out=[];
out = 'Goodness of fit: r2(Vel/WithoutDC2) | r2(Vel/WithDC2) '; 
out = strvcat(out, sprintf('--------------------------------------------------------------------------------------------------------'));    
clear OutputValue;OutputValue=[r_squared_Vel_NDC2 r_squared_Vel];
out=strvcat(out, sprintf('                   %6.3f           | %6.3f ', OutputValue));
figure(FigureIndex);axes('position',[0.26 0.65 0.25 0.25]); set(gca,'box','off','visible','off');
text(-0.08,1,out,'fontsize',8,'fontname','courier','horizontalalignment','left','verticalalignment','top');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% order in which the directions are plotted
plot_col = [1 1 1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 6 6 6 7 7 7 8 8 8];
plot_row = [5 4 3 2 1 4 3 2 4 3 2 4 3 2 4 3 2 4 3 2 4 3 2 4 3 2];

x_start = [0, 0];
x_stop =  [2, 2];
y_marker=[0,1.1*max(max(spacetime_data))];

xscale = [0.8 0.77 0.65 0.52 0.5 0.52 0.65 0.77];   
yscale = [0.35 0.5 0.65 0.5 0.35 0.2 0.05 0.2];

%get the curve fitting results for each direction
clear model_Vel; model_Vel=zeros(length(Azi_temp),length(x_time));
model_Vel(:,1:StartIndex-1) = NaN; model_Vel(:,StartIndex:length(x_time))=spacefit_Vel(1:8,:);

clear model_Vel_NDC2; model_Vel_NDC2=zeros(length(Azi_temp),length(x_time));
model_Vel_NDC2(:,1:StartIndex-1) = NaN; model_Vel_NDC2(:,StartIndex:length(x_time))=spacefit_Vel_NDC2(1:8,:);

%Plot the PSTH, superimposed with curve fitting
for i=1:length(Azi_temp)
    i;
    axes('position',[xscale(i) yscale(i) 0.16 0.10]);
    bar(x_time,PSTH_tempPlane(i,:),1.0);hold on;%bar(x_time, count_y{i,j,k});    hold on;
    
    plot(x_time,model_Vel(i,:),'r','LineWidth',2); hold on;   
    plot(x_time,model_Vel_NDC2(i,:),'g','LineWidth',2); hold on;   
    
    plot( x_start, y_marker, 'k-','LineWidth',2);hold on;
    plot( x_stop,  y_marker, 'k-','LineWidth',2);hold on;
    xlim([-0.5,2.5]);    ylim([0,1.2*max(max(spacetime_data))]);
    
    set(gca, 'xtick', [] );set(gca, 'XTickMode','manual'); 
    set(gca, 'xtick',[-0.5:0.5:2.5]);
    set( gca, 'xticklabel', '-0.5|0|0.5|1|1.5|2|2.5');
    set( gca, 'yticklabel', ' ' );
    title({['azi=' num2str(Azi_temp(i)) '; Ele=' num2str(Ele_temp(i))]}, 'FontSize', 8);
end  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output some data for further analysis
CurveFit.spacetime_data=spacetime_data;
CurveFit.PeakValue=PeakValue_tempPlane;
CurveFit.TroughValue=TroughValue_tempPlane;
CurveFit.p_peak=p_peak_tempPlane;
CurveFit.p_trough=p_trough_tempPlane;
CurveFit.p_anova_2way=p_anova_2way;

CurveFit.spacefit_Vel=spacefit_Vel;
CurveFit.vect_Vel=vect_Vel;
CurveFit.CI_Vel=CI_Vel;
CurveFit.cor_Vel=cor_Vel;
CurveFit.err_Vel=err_Vel;
CurveFit.r_squared_Vel=r_squared_Vel;

CurveFit.spacefit_Vel_NDC2=spacefit_Vel_NDC2;
CurveFit.vect_Vel_NDC2=vect_Vel_NDC2;
CurveFit.CI_Vel_NDC2=CI_Vel_NDC2;
CurveFit.cor_Vel_NDC2=cor_Vel_NDC2;
CurveFit.err_Vel_NDC2=err_Vel_NDC2;
CurveFit.r_squared_Vel_NDC2=r_squared_Vel_NDC2;


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


