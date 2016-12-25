function RF_SUMU(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

ProtocolDefs;
Path_Defs;
SpikeChan=1

RFContourFigurePlot(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Directory ='C:\Aihua\z_TempOutputs\figures';
FileName1=[Directory,'\',FILE(1:end-4),'_SU_01_Direction_RF.fig'];figure(2); saveas(gcf,FileName1,'fig');
FileName2=[Directory,'\',FILE(1:end-4),'_SU_02_Direction_RF_forward.fig'];figure(3); saveas(gcf,FileName2,'fig');
FileName3=[Directory,'\',FILE(1:end-4),'_SU_03_Wrapped Gaussian Fitting.fig'];figure(4); saveas(gcf,FileName3,'fig');
FileName4=[Directory,'\',FILE(1:end-4),'_SU_04_Goodness of Fitting.fig'];figure(5); saveas(gcf,FileName4,'fig');
FileName5=[Directory,'\',FILE(1:end-4),'_SU_05_CirVar.fig'];figure(6); saveas(gcf,FileName5,'fig');
FileName6=[Directory,'\',FILE(1:end-4),'_SU_06_rVar.fig'];figure(7);saveas(gcf,FileName6,'fig');
FileName7=[Directory,'\',FILE(1:end-4),'_SU_07_Spatial_RF.fig'];figure(8);saveas(gcf,FileName7,'fig');

% SpikeChan=4
% RFContourFigurePlot(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Directory ='C:\Aihua\z_TempOutputs\figures';
% FileName1=[Directory,'\',FILE(1:end-4),'_MU_01_Direction_RF.fig'];figure(2); saveas(gcf,FileName1,'fig');
% FileName2=[Directory,'\',FILE(1:end-4),'_MU_02_Direction_RF_forward.fig'];figure(3); saveas(gcf,FileName2,'fig');
% FileName3=[Directory,'\',FILE(1:end-4),'_MU_03_Wrapped Gaussian Fitting.fig'];figure(4); saveas(gcf,FileName3,'fig');
% FileName4=[Directory,'\',FILE(1:end-4),'_MU_04_Goodness of Fitting.fig'];figure(5); saveas(gcf,FileName4,'fig');
% FileName5=[Directory,'\',FILE(1:end-4),'_MU_05_CirVar.fig'];figure(6); saveas(gcf,FileName5,'fig');
% FileName6=[Directory,'\',FILE(1:end-4),'_MU_06_rVar.fig'];figure(7);saveas(gcf,FileName6,'fig');
% FileName7=[Directory,'\',FILE(1:end-4),'_MU_07_Spatial_RF.fig'];figure(8);saveas(gcf,FileName7,'fig');



return;