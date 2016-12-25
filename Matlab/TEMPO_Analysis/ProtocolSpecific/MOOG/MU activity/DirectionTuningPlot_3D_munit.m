% DirectionTuningPlot_3D.m -- Plots response as a function of azimuth and elevation for MOOG 3D tuning expt
%--	YONG, 6/27/03
%-----------------------------------------------------------------------------------------------------------------------
function DirectionTuningPlot_3D(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

Path_Defs;
ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP
%---------------------------------------------------------------------------------------
% PeakValue_SU=max(data.spike_rates(1, :));
% PeakValue_MU=max(data.spike_rates(4, :));
% %Also, write out some summary data to a cumulative summary file
% buff = sprintf('%s\t %4.2f\t %4.3f\t %4.3f\t',FILE,PeakValue_SU,PeakValue_MU);
% %outfile = [BASE_PATH 'ProtocolSpecific\MOOG\MU activity\Corr.dat'];
% outfile = ['Z:\Users\Aihua\Analysis\Analysis_SU_MU\data\PeakValue.dat'];
% printflag = 0;
% if (exist(outfile, 'file') == 0)    %file does not yet exist
%     printflag = 1;
% end
% fid = fopen(outfile, 'a');
% if (printflag)
%     fprintf(fid, 'FILE\t         PeakValue_SU\t PeakValue_MU\t');
%     fprintf(fid, '\r\n');
% end
% fprintf(fid, '%s', buff);
% fprintf(fid, '\r\n');
% fclose(fid);



figure(2);legend('Before','SU Removed'); Title=[FILE ' // Cross correlation between SU and MU'];title(Title);

Ori=data.corr{1};
Removed=data.corr{2};

X1=round(length(Ori)/2)-1:round(length(Ori)/2)+1;
X2=1:3;
X3=length(Ori)-2:length(Ori);

Before_A1 = trapz(X1,Ori(X1))
Before_A2 = trapz(X2,Ori(X2))+trapz(X3,Ori(X3))
Before=Before_A1/Before_A2

After_A1 = trapz(X1,Removed(X1))
After_A2 = trapz(X2,Removed(X2))+trapz(X3,Removed(X3))
After=After_A1/After_A2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save_name = 'CorrOri'; % name of the  mat file with all the Original CrossCorrelation Histrogram 
dum_name = FILE(1:end-4);
spacer = ' ';
eval([dum_name '= Ori;'])
if (exist([save_name '.mat'], 'file') == 0)    %file does not yet exist 
    eval(['save Z:\Users\Aihua\Cluster_Analysis\Analysis_SU_MU\data\' save_name spacer dum_name]);   
    %Z:\Users\Aihua\Cluster_Analysis\Analysis_SU_MU\data
else 
    eval(['save Z:\Users\Aihua\Cluster_Analysis\Analysis_SU_MU\data\' save_name spacer dum_name ' -APPEND']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save_name = 'CorrAfter'; % name of the  mat file with all the CrossCorrelation Histrogram after SU Removed
dum_name = FILE(1:end-4);
spacer = ' ';
eval([dum_name '= Removed;'])
if (exist([save_name '.mat'], 'file') == 0)    %file does not yet exist  
    eval(['save Z:\Users\Aihua\Cluster_Analysis\Analysis_SU_MU\data\' save_name spacer dum_name]);
else
    eval(['save Z:\Users\Aihua\Cluster_Analysis\Analysis_SU_MU\data\' save_name spacer dum_name ' -APPEND']);
end


%---------------------------------------------------------------------------------------
%Also, write out some summary data to a cumulative summary file
buff = sprintf('%s\t %4.2f\t %4.3f\t %4.3f\t %4.3f\t %4.3f\t %4.3f\t',FILE,Before_A1,Before_A2,After_A1,After_A2,Before, After );
%outfile = [BASE_PATH 'ProtocolSpecific\MOOG\MU activity\Corr.dat'];
outfile = ['Z:\Users\Aihua\Cluster_Analysis\Analysis_SU_MU\data\Corr.dat'];
printflag = 0;
if (exist(outfile, 'file') == 0)    %file does not yet exist
    printflag = 1;
end
fid = fopen(outfile, 'a');
if (printflag)
    fprintf(fid, 'FILE\t         Before_A1\t Before_A2\t After_A1\t After_A2\t Before\t After\t');
    fprintf(fid, '\r\n');
end
fprintf(fid, '%s', buff);
fprintf(fid, '\r\n');
fclose(fid);


Directory='Z:\Users\Aihua\Cluster_Analysis\Analysis_SU_MU\figure';
FileName0=[Directory,'\',FILE(1:end-4),'_corr.fig'];figure(2); saveas(gcf,FileName0,'fig');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot SU activity
FigureIndex=2;
TuningPlotNew(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE,FigureIndex);

%plot MU activity
FigureIndex=3;
SpikeChan=4;
TuningPlotNew(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE,FigureIndex);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------------------------------------------------------------------
%save the figures
Directory='Z:\Users\Aihua\Cluster_Analysis\Analysis_SU_MU\figure';
FileName1=[Directory,'\',FILE(1:end-4),'_SU.fig'];figure(2); saveas(gcf,FileName1,'fig');
FileName2=[Directory,'\',FILE(1:end-4),'_MU.fig'];figure(3); saveas(gcf,FileName2,'fig');

return;