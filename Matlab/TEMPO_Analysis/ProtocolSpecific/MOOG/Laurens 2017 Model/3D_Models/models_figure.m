% plot figures of 3D models for PSTH
% red: velocity only (VO)
% blue: acceleration only (AO)
% green: VA
% magenta: VJ
% cyan: AJ
% black: VAJ
% white: PVAJ
% 100# PSTH model figures
% 110# contour model figures
% LBY 20170328
% LBY 20170603
% LBY 20171210 add contour figures

global PSTH3Dmodel PSTH;

nBins = stimOffBin - stimOnBin +1;
%
% k=1,2,3
% j= -90,-45,0,45,90 (up->down)
% i=0 45 90 135 180 225 270 315

stimType{1}='Vestibular';
stimType{2}='Visual';
stimType{3}='Combined';

% 270-225-180-135-90-45-0-315-270 for figures
iAzi = [7 6 5 4 3 2 1 8 7];
xAzi = 1:9;
yEle = [-1,-0.707,0,0.707,1];
% initialize default properties
set(0,'defaultaxesfontsize',24);
colorDefsLBY;

% the lines markers
markers = {
    'aMax',stimOnBin+aMax/timeStep-1,colorDBlue;
    'aMin',stimOnBin+aMin/timeStep-1,colorDBlue;
    'vMax',stimOnBin+(aMin+aMax)/2/timeStep-1,colorDRed;
    };

unique_elevation = [-90,-45,0,45,90];
unique_azimuth = [0 45 90 135 180 225 270 315];

sprintf('Plot model figures...')

% ------ fig.100 plot 3D models (PSTH) ------%

for m_inx = 1:length(models)
    figure(100+m_inx);
    set(gcf,'pos',[0 0 1900 1000]);
    clf;
    [~,h_subplot] = tight_subplot(5,9,0.04,0.15);
    
    for j = 2:length(unique_elevation)-1
        for i = 1:length(unique_azimuth)+1
            axes(h_subplot(i+(j-1)*9));hold on;
            bar(squeeze(PSTH_data(j,iAzi(i),:)),'facecolor',colorLGray,'edgecolor',colorLGray);
            for n = 1:size(markers,1)
                plot([markers{n,2} markers{n,2}], [0,max(PSTH.maxSpkRealBinMean(stimTypeInx),PSTH.maxSpkSponBinMean)], '--','color',markers{n,3},'linewidth',2);
                hold on;
            end
            eval(['h = plot(1:nBins,squeeze(PSTH3Dmodel{',num2str(stimTypeInx),'}.modelFitRespon_',models{m_inx},'(',num2str(j),',',num2str(iAzi(i)),',:)));']);
            %             set(h,'linestyle','-','marker','.','color',models_color{m_inx},'markersize',8);
            set(h,'linestyle','-','linewidth',2,'color',models_color{m_inx});
            %                 keyboard;
            set(gca,'ylim',[0 max(max(PSTH_data(:)),max(temp_spon))],'xlim',[1 nBins]);
            axis on;
            SetFigure(15);
            set(gca,'xtick',[],'xticklabel',[]);
            
        end
    end
    
    % 2 extra conditions
    axes(h_subplot(5+(1-1)*9));hold on;
    bar(squeeze(PSTH_data(1,iAzi(5),:)),'facecolor',colorLGray,'edgecolor',colorLGray);
    for n = 1:size(markers,1)
        plot([markers{n,2} markers{n,2}], [0,max(PSTH.maxSpkRealBinMean(stimTypeInx),PSTH.maxSpkSponBinMean)], '--','color',markers{n,3},'linewidth',2);
        hold on;
    end
    eval(['h = plot(1:nBins,squeeze(PSTH3Dmodel{',num2str(stimTypeInx),'}.modelFitRespon_',models{m_inx},'(',num2str(1),',',num2str(iAzi(5)),',:)));']);
    %     set(h,'linestyle','none','marker','.','color',models_color{m_inx},'markersize',8);
    set(h,'linestyle','-','linewidth',2,'color',models_color{m_inx});
    set(gca,'ylim',[0 max(max(PSTH_data(:)),max(temp_spon))],'xlim',[1 nBins]);
    axis on;
    SetFigure(15);
    set(gca,'xtick',[],'xticklabel',[]);
    
    
    axes(h_subplot(5+(5-1)*9));hold on;
    bar(squeeze(PSTH_data(5,iAzi(5),:)),'facecolor',colorLGray,'edgecolor',colorLGray);
    for n = 1:size(markers,1)
        plot([markers{n,2} markers{n,2}], [0,max(PSTH.maxSpkRealBinMean(stimTypeInx),PSTH.maxSpkSponBinMean)], '--','color',markers{n,3},'linewidth',2);
        hold on;
    end
    eval(['h = plot(1:nBins,squeeze(PSTH3Dmodel{',num2str(stimTypeInx),'}.modelFitRespon_',models{m_inx},'(',num2str(5),',',num2str(iAzi(5)),',:)));']);
    %     set(h,'linestyle','none','marker','.','color',models_color{m_inx},'markersize',8);
    set(h,'linestyle','-','linewidth',2,'color',models_color{m_inx});
    set(gca,'ylim',[0 max(max(PSTH_data(:)),max(temp_spon))],'xlim',[1 nBins]);
    axis on;
    SetFigure(15);
    set(gca,'xtick',[],'xticklabel',[]);
    
    
    % sponteneous
    axes(h_subplot(1+(1-1)*9));hold on;
    bar(temp_spon,'facecolor',colorLGray,'edgecolor',colorLGray);
    %     for n = 1:size(markers,1)
    %         plot([markers{n,2} markers{n,2}], [0,max(PSTH.maxSpkRealBinMean(stimTypeInx),PSTH.maxSpkSponBinMean)], '--','color',markers{n,3},'linewidth',0.5);
    %         hold on;
    %     end
    set(gca,'ylim',[0 max(max(PSTH_data(:)),max(temp_spon))],'xlim',[1 nBins]);
    axis on;
    SetFigure(15);
    set(gca,'xtick',[],'xticklabel',[]);
    
    
    % text on the figure
    axes('unit','pixels','pos',[60 800 1800 180]);
    xlim([0,100]);
    ylim([0,20]);
    switch Protocol
        case 100
            text(40,14,'PSTHs for each direction(T)','fontsize',20);
        case 112
            text(40,14,'PSTHs for each direction(R)','fontsize',20);
    end
    text(1,10,'Spontaneous','fontsize',15);
    FileNameTemp = num2str(FILE);
    %     FileNameTemp =  FileNameTemp(1:end-4);
    str = FileNameTemp;
    str1 = [FileNameTemp,'  ',models{m_inx},' model_','    ',stimType{stimTypeInx}];
    text(60,6,str1,'fontsize',18);
    eval(['R_squared = ','PSTH3Dmodel{',num2str(stimTypeInx),'}.RSquared_',models{m_inx},';']);
    text(60,0,['R squared = ',num2str(R_squared)],'fontsize',18);
    axis off;
    
    axes('unit','pixels','pos',[60 50 1800 100]);
    xlim([0,100]);
    ylim([0,10]);
    %     text(0,0,['Max spk mean FR(Hz): ',num2str(maxSpkRealAll(k))],'fontsize',15);
    text(0,10,['Spon max mean FR(Hz): ',num2str(max(temp_spon))],'fontsize',15);
    text(0,6,['Spon mean FR(Hz): ',num2str(meanSpon)],'fontsize',15);
    
    %--- this is for annotation - the direction ----%
    %     text(3,0,'\downarrow ','fontsize',30);
    %     text(26,0,'\leftarrow ','fontsize',30);
    %     text(49,0,'\uparrow ','fontsize',30);
    %     text(71,0,'\rightarrow ','fontsize',30);
    %     text(93,0,'\downarrow ','fontsize',30);
    %     text(55,76,'\uparrow ','fontsize',30);
    %     text(55,13,'\downarrow ','fontsize',30);
    axis off;
    %--- this is for annotation - the direction ----%
    
    
    % to save the figures
    str3 = [str, '_',models{m_inx},'_model_',stimType{stimTypeInx}];
    set(gcf,'paperpositionmode','auto');
    if Protocol == 100
        ss = [str3, '_T'];
%         saveas(100+m_inx,['Z:\LBY\Recording data\Qiaoqiao\3D_Tuning_models\Translation\' ss], 'emf');
saveas(100+m_inx,['Z:\LBY\Recording data\Zebulon\MSTd\Models\' ss], 'emf'); % for MSTd 
        
    elseif Protocol == 112
        ss = [str3, '_R'];
%         saveas(100+m_inx,['Z:\LBY\Recording data\Qiaoqiao\3D_Tuning_models\Rotation\' ss], 'emf');
    end
end


% % ------ fig.110 plot 3D models (Contour) ------%
%
%
% for m_inx = 1:length(models)
%     figure(110+m_inx);
%     set(gcf,'pos',[60 200 1000 600]);
%     clf;
%
%     axes('unit','pixels','pos',[200 100 650 300]);
%     %     contourf(xAzi,yEle,spk_data_count_mean_rate_trans{k},'linecolor','w','linestyle','none');
%     eval(['h = contourf(xAzi,yEle,PSTH3Dmodel{',num2str(stimTypeInx),'}.modelFitResponMeanTrans_',models{m_inx},');']);
%     colorbar;
%     set(gca, 'ydir' , 'reverse'); % so that up is up, down is down
%     set(gca, 'xtick', [] );
%     set(gca, 'ytick', [] );
%     box off;
%
%     % text on the figure
%     axes('unit','pixels','pos',[60 470 1000 100]);
%     xlim([0,100]);
%     ylim([0,20]);
%     switch Protocol
%         case 100
%             text(15,14,'Color contour maps - model fitting(T)','fontsize',20);
%         case 112
%             text(15,14,'Color contour maps - model fitting(R)','fontsize',20);
%     end
%     FileNameTemp = num2str(FILE);
%     %     FileNameTemp =  FileNameTemp(1:end-4);
%     str = FileNameTemp;
%     str1 = [FileNameTemp,'  ',models{m_inx},' model_','    ',stimType{stimTypeInx}];
%     text(40,4,str1,'fontsize',18);
%     eval(['R_squared = ','PSTH3Dmodel{',num2str(stimTypeInx),'}.RSquared_',models{m_inx},';']);
%     text(40,-2,['R squared = ',num2str(R_squared)],'fontsize',18);
%     axis off;
%
%     axes('unit','pixels','pos',[60 50 1000 100]);
%     xlim([0,100]);
%     ylim([0,10]);
%     text(0,3,['Spon max mean FR(Hz): ',num2str(max(temp_spon))],'fontsize',15);
%     text(0,-1,['Spon mean FR(Hz): ',num2str(meanSpon)],'fontsize',15);
%     axis off;
%
%     % to save the figures
%     str3 = [str, '_',models{m_inx},'_model_',stimType{stimTypeInx}];
%     set(gcf,'paperpositionmode','auto');
%     if Protocol == 100
%         ss = [str3, '_Contour_T'];
%         saveas(110+m_inx,['Z:\LBY\Recording data\Qiaoqiao\3D_Tuning_models\Translation\' ss], 'emf');
%     elseif Protocol == 112
%         ss = [str3, '_Contour_R'];
%         saveas(110+m_inx,['Z:\LBY\Recording data\Qiaoqiao\3D_Tuning_models\Rotation\' ss], 'emf');
%     end
%
% end

% % ------ fig.120 plot V,A,J direction tuning for VAJ models (Contour) ------%
%
%
% components = {'V','A','J'};
% figure(120);clf;set(gcf,'pos',[20 200 1800 200]);
% for inx = 1:length(components)
%
%     subplot(1,length(components),inx)
%     eval(['h = contourf(xAzi,yEle,PSTH3Dmodel{',num2str(stimTypeInx),'}.modelFitMeanTrans_VAJ_',components{inx},');']);
%     colorbar;
% %     set(gca, 'linecolor' , 'w');
%     set(gca, 'ydir' , 'reverse'); % so that up is up, down is down
%     set(gca, 'xtick', [] );
%     set(gca, 'ytick', [] );
%     box off;
% end
% % text on the figure
% axes('unit','pixels','pos',[60 270 1800 100]);
% xlim([0,100]);
% ylim([0,20]);
% switch Protocol
%     case 100
%         text(15,14,'Color contour maps - model fitting(T)','fontsize',20);
%     case 112
%         text(15,14,'Color contour maps - model fitting(R)','fontsize',20);
% end
% FileNameTemp = num2str(FILE);
% %     FileNameTemp =  FileNameTemp(1:end-4);
% str = FileNameTemp;
% str1 = [FileNameTemp,'  ',components{inx},' model_','    ',stimType{stimTypeInx}];
% text(40,4,str1,'fontsize',18);
% eval(['R_squared = ','PSTH3Dmodel{',num2str(stimTypeInx),'}.RSquared_VAJ;']);
% text(40,-2,['R squared = ',num2str(R_squared)],'fontsize',18);
% axis off;
%
% % axes('unit','pixels','pos',[60 50 1000 100]);
% % xlim([0,100]);
% % ylim([0,10]);
% % text(0,3,['Spon max mean FR(Hz): ',num2str(max(temp_spon))],'fontsize',15);
% % text(0,-1,['Spon mean FR(Hz): ',num2str(meanSpon)],'fontsize',15);
% % axis off;
%
% % to save the figures
% str3 = [str, '_VAJmodel_',stimType{stimTypeInx}];
% set(gcf,'paperpositionmode','auto');
% if Protocol == 100
%     ss = [str3, '_Contour_T'];
%     saveas(120,['Z:\LBY\Recording data\Qiaoqiao\3D_Tuning_models\Translation\' ss], 'emf');
% elseif Protocol == 112
%     ss = [str3, '_Contour_R'];
%     saveas(120,['Z:\LBY\Recording data\Qiaoqiao\3D_Tuning_models\Rotation\' ss], 'emf');
% end
%
% % ------ fig.130 plot V,A,J direction tuning for VAJ models at corresponding peak time(Contour) ------%
%
%
% components = {'V','A','J'};
% figure(130);clf;set(gcf,'pos',[20 200 1800 200]);
% for inx = 1:length(components)
%
%     subplot(1,length(components),inx)
%     eval(['h = contourf(xAzi,yEle,PSTH3Dmodel{',num2str(stimTypeInx),'}.modelFitMeanTrans_VAJ_',components{inx},'_peakT);']);
%     colorbar;
% %     set(gca, 'linecolor' , 'w');
%     set(gca, 'ydir' , 'reverse'); % so that up is up, down is down
%     set(gca, 'xtick', [] );
%     set(gca, 'ytick', [] );
%     box off;
% end
% % text on the figure
% axes('unit','pixels','pos',[60 270 1800 100]);
% xlim([0,100]);
% ylim([0,20]);
% switch Protocol
%     case 100
%         text(15,14,'Color contour maps - model fitting(T)','fontsize',20);
%     case 112
%         text(15,14,'Color contour maps - model fitting(R)','fontsize',20);
% end
% FileNameTemp = num2str(FILE);
% %     FileNameTemp =  FileNameTemp(1:end-4);
% str = FileNameTemp;
% str1 = [FileNameTemp,'  VAJ model_','    ',stimType{stimTypeInx}];
% text(40,4,str1,'fontsize',18);
% eval(['R_squared = ','PSTH3Dmodel{',num2str(stimTypeInx),'}.RSquared_VAJ;']);
% text(40,-2,['R squared = ',num2str(R_squared)],'fontsize',18);
% axis off;
%
% % axes('unit','pixels','pos',[60 50 1000 100]);
% % xlim([0,100]);
% % ylim([0,10]);
% % text(0,3,['Spon max mean FR(Hz): ',num2str(max(temp_spon))],'fontsize',15);
% % text(0,-1,['Spon mean FR(Hz): ',num2str(meanSpon)],'fontsize',15);
% % axis off;
%
% % to save the figures
% str3 = [str, '_VAJmodel_peakT_',stimType{stimTypeInx}];
% set(gcf,'paperpositionmode','auto');
% if Protocol == 100
%     ss = [str3, '_Contour_T'];
%     saveas(130,['Z:\LBY\Recording data\Qiaoqiao\3D_Tuning_models\Translation\' ss], 'emf');
% elseif Protocol == 112
%     ss = [str3, '_Contour_R'];
%     saveas(130,['Z:\LBY\Recording data\Qiaoqiao\3D_Tuning_models\Rotation\' ss], 'emf');
% end

% % ------ fig.140 plot V,A spatial tuning for VA models (Contour, use fitted data) ------%
% 
% 
% components = {'V','A'};
% figure(140);clf;set(gcf,'pos',[120 200 1400 600]);
% [~,h_subplot] = tight_subplot(1,2,0.1,0.3);
% for inx = 1:length(components)
%     
%     axes(h_subplot(inx));
%     eval(['h = contourf(xAzi,yEle,PSTH3Dmodel{',num2str(stimTypeInx),'}.modelFitTrans_spatial_VA.',components{inx},');']);
%     a = cell2mat(struct2cell(PSTH3Dmodel{stimTypeInx}.modelFitTrans_spatial_VA));
%     caxis([min(min(a)) max(max(a))]);
%     colorbar;
%     %     set(gca, 'linecolor' , 'w');
%     set(gca, 'ydir' , 'reverse'); % so that up is up, down is down
%     set(gca, 'xtick', [] );
%     set(gca, 'ytick', [] );
%     box off;
%     title(['Spatial tuning - ',components{inx}]);
% end
% % text on the figure
% axes('unit','pixels','pos',[60 470 1400 100]);
% xlim([0,100]);
% ylim([0,20]);
% switch Protocol
%     case 100
%         text(11,16,'Color contour maps - model fitting(T)','fontsize',20);
%     case 112
%         text(11,16,'Color contour maps - model fitting(R)','fontsize',20);
% end
% FileNameTemp = num2str(FILE);
% %     FileNameTemp =  FileNameTemp(1:end-4);
% str = FileNameTemp;
% str1 = [FileNameTemp,'  VA model_','    ',stimType{stimTypeInx}];
% text(50,14,str1,'fontsize',18);
% eval(['R_squared = ','PSTH3Dmodel{',num2str(stimTypeInx),'}.RSquared_VA;']);
% text(60,7,['R squared = ',num2str(R_squared)],'fontsize',18);
% axis off;
% 
% % axes('unit','pixels','pos',[60 50 1000 100]);
% % xlim([0,100]);
% % ylim([0,10]);
% % text(0,3,['Spon max mean FR(Hz): ',num2str(max(temp_spon))],'fontsize',15);
% % text(0,-1,['Spon mean FR(Hz): ',num2str(meanSpon)],'fontsize',15);
% % axis off;
% 
% % to save the figures
% str3 = [str, '_VAmodel_',stimType{stimTypeInx}];
% set(gcf,'paperpositionmode','auto');
% if Protocol == 100
%     ss = [str3, '_Contour_T'];
%     saveas(140,['Z:\LBY\Recording data\Qiaoqiao\3D_Tuning_models\Translation\' ss], 'emf');
% elseif Protocol == 112
%     ss = [str3, '_Contour_R'];
%     saveas(140,['Z:\LBY\Recording data\Qiaoqiao\3D_Tuning_models\Rotation\' ss], 'emf');
% end
% 
% ------ fig.150 plot V,A for VA model (PSTH) ------%
components = {'V','A'};



for m_inx = 1:length(components)
%     figure(150+m_inx);clf;set(gcf,'pos',[20+(m_inx-1)*910 200 900 500]);
figure(150+m_inx);clf;set(gcf,'pos',[0 0 1900 1000]);
    [~,h_subplot] = tight_subplot(5,9,0.04,0.15);
    
    for j = 2:length(unique_elevation)-1
        for i = 1:length(unique_azimuth)+1
            axes(h_subplot(i+(j-1)*9));hold on;
            bar(squeeze(PSTH_data(j,iAzi(i),:)),'facecolor',colorLGray,'edgecolor',colorLGray);
            for n = 1:size(markers,1)
                plot([markers{n,2} markers{n,2}], [0,max(PSTH.maxSpkRealBinMean(stimTypeInx),PSTH.maxSpkSponBinMean)], '--','color',markers{n,3},'linewidth',2);
                hold on;
            end
            eval(['h = plot(1:nBins,squeeze(PSTH3Dmodel{',num2str(stimTypeInx),'}.modelFitMeanTrans_VA.',components{m_inx},'(',num2str(j),',',num2str(iAzi(i)),',:)));']);
            %             set(h,'linestyle','-','marker','.','color',models_color{m_inx},'markersize',8);
            set(h,'linestyle','-','linewidth',2,'color','k');
            %                 keyboard;
            set(gca,'ylim',[0 max(max(PSTH_data(:)),max(temp_spon))],'xlim',[1 nBins]);
            axis on;
            
            set(gca,'xtick',[],'xticklabel',[]);
            
        end
    end
    
    % 2 extra conditions
    axes(h_subplot(5+(1-1)*9));hold on;
    bar(squeeze(PSTH_data(1,iAzi(5),:)),'facecolor',colorLGray,'edgecolor',colorLGray);
    for n = 1:size(markers,1)
        plot([markers{n,2} markers{n,2}], [0,max(PSTH.maxSpkRealBinMean(stimTypeInx),PSTH.maxSpkSponBinMean)], '--','color',markers{n,3},'linewidth',2);
        hold on;
    end
    eval(['h = plot(1:nBins,squeeze(PSTH3Dmodel{',num2str(stimTypeInx),'}.modelFitMeanTrans_VA.',components{m_inx},'(',num2str(1),',',num2str(iAzi(5)),',:)));']);
    %     set(h,'linestyle','none','marker','.','color',models_color{m_inx},'markersize',8);
    set(h,'linestyle','-','linewidth',2,'color','k');
    set(gca,'ylim',[0 max(max(PSTH_data(:)),max(temp_spon))],'xlim',[1 nBins]);
    axis on;
%     SetFigure(15);
    set(gca,'xtick',[],'xticklabel',[]);
    
    
    axes(h_subplot(5+(5-1)*9));hold on;
    bar(squeeze(PSTH_data(5,iAzi(5),:)),'facecolor',colorLGray,'edgecolor',colorLGray);
    for n = 1:size(markers,1)
        plot([markers{n,2} markers{n,2}], [0,max(PSTH.maxSpkRealBinMean(stimTypeInx),PSTH.maxSpkSponBinMean)], '--','color',markers{n,3},'linewidth',2);
        hold on;
    end
    eval(['h = plot(1:nBins,squeeze(PSTH3Dmodel{',num2str(stimTypeInx),'}.modelFitMeanTrans_VA.',components{m_inx},'(',num2str(5),',',num2str(iAzi(5)),',:)));']);
    %     set(h,'linestyle','none','marker','.','color',models_color{m_inx},'markersize',8);
    set(h,'linestyle','-','linewidth',2,'color','k');
    set(gca,'ylim',[0 max(max(PSTH_data(:)),max(temp_spon))],'xlim',[1 nBins]);
    axis on;
%     SetFigure(15);
    set(gca,'xtick',[],'xticklabel',[]);
    
    
    % sponteneous
    axes(h_subplot(1+(1-1)*9));hold on;
    bar(temp_spon,'facecolor',colorLGray,'edgecolor',colorLGray);
    %     for n = 1:size(markers,1)
    %         plot([markers{n,2} markers{n,2}], [0,max(PSTH.maxSpkRealBinMean(stimTypeInx),PSTH.maxSpkSponBinMean)], '--','color',markers{n,3},'linewidth',0.5);
    %         hold on;
    %     end
    set(gca,'ylim',[0 max(max(PSTH_data(:)),max(temp_spon))],'xlim',[1 nBins]);
    axis on;
%     SetFigure(15);
    set(gca,'xtick',[],'xticklabel',[]);
    
    
    % text on the figure
    axes('unit','pixels','pos',[60 800 1800 100]);
    xlim([0,100]);
    ylim([0,20]);
    switch Protocol
        case 100
            text(30,20,'PSTH (T)','fontsize',20);
        case 112
            text(30,20,'PSTH (R)','fontsize',20);
    end
    text(10,2,'Spontaneous','fontsize',15);
    FileNameTemp = num2str(FILE);
    %     FileNameTemp =  FileNameTemp(1:end-4);
    str = FileNameTemp;
    str1 = [FileNameTemp,'  VA model ',components{m_inx},'  ',stimType{stimTypeInx}];
    text(40,20,str1,'fontsize',18);
    eval(['R_squared = ','PSTH3Dmodel{',num2str(stimTypeInx),'}.RSquared_VA;']);
    text(70,10,['R squared = ',num2str(R_squared)],'fontsize',18);
    axis off;
    
    %     axes('unit','pixels','pos',[60 50 1800 100]);
    %     xlim([0,100]);
    %     ylim([0,10]);
    %     %     text(0,0,['Max spk mean FR(Hz): ',num2str(maxSpkRealAll(k))],'fontsize',15);
    %     text(0,10,['Spon max mean FR(Hz): ',num2str(max(temp_spon))],'fontsize',15);
    %     text(0,6,['Spon mean FR(Hz): ',num2str(meanSpon)],'fontsize',15);
    
    %--- this is for annotation - the direction ----%
    %     text(3,0,'\downarrow ','fontsize',30);
    %     text(26,0,'\leftarrow ','fontsize',30);
    %     text(49,0,'\uparrow ','fontsize',30);
    %     text(71,0,'\rightarrow ','fontsize',30);
    %     text(93,0,'\downarrow ','fontsize',30);
    %     text(55,76,'\uparrow ','fontsize',30);
    %     text(55,13,'\downarrow ','fontsize',30);
    axis off;
    %--- this is for annotation - the direction ----%
    
    SetFigure(15);
    % to save the figures
    str3 = [str, '_VA_model_',components{m_inx},'_',stimType{stimTypeInx}];
    set(gcf,'paperpositionmode','auto');
    if Protocol == 100
        ss = [str3, '_T'];
        saveas(150+m_inx,['Z:\LBY\Recording data\Qiaoqiao\3D_Tuning_models\Translation\' ss], 'emf');
    elseif Protocol == 112
        ss = [str3, '_R'];
        saveas(150+m_inx,['Z:\LBY\Recording data\Qiaoqiao\3D_Tuning_models\Rotation\' ss], 'emf');
    end
end

toc;
