%% plot figures - PSTH
%
global PSTH;

% ------ fig.20 plot PSTH for each trial and raster plot across directions ------%

for k = 1:length(unique_stimType)
    figure(20+k);
    set(gcf,'pos',[0 0 1900 1000]);
    clf;
    for j = 2:length(unique_elevation)-1
        for i = 1:length(unique_azimuth)+1
            h1 = axes('unit','pixels','pos',[30+200*(i-1) 20+180*(5-j) 180 80]);
            for n = 1:size(markers,1)
                plot(h1,[markers{n,3} markers{n,3}], [0,max(PSTH.maxSpkRealBinMean(k),PSTH.maxSpkSponBinMean)], '--','color',markers{n,4},'linewidth',0.5);
                hold on;
            end
            for r = 1:size((PSTH.spk_data_bin_rate{k,j,iAzi(i)}),2)
                color = ['color',num2str(r)];
                plot(h1,PSTH.spk_data_bin_rate{k,j,iAzi(i)}(:,r)','color',eval(color));
                axis off;
                hold on;
                h2 = axes('unit','pixels','pos',[30+200*(i-1) 20+90+180*(5-j)+8*(r-1) 180 8]);
                %                 tt = stimOnT(1)-tOffset1:stimOffT(1)+tOffset2; yy = ones(tOffset1+tOffset2+unique_duration(1,1));
                %                 tt = tt(logical(spk_data{k,j,iAzi(i)}(stimOnT(1)-tOffset1:stimOffT(1)+tOffset2,r)));
                %                 yy = yy(logical(spk_data{k,j,iAzi(i)}(stimOnT(1)-tOffset1:stimOffT(1)+tOffset2,r)));
                %                 yy = repmat(yy,1,2);
                %                 plot(h2,tt,yy,'.k');
                stem(h2, stimOnT(1)-tOffset1:stimOffT(1)+tOffset2,spk_data{k,j,iAzi(i)}(stimOnT(1)-tOffset1:stimOffT(1)+tOffset2,r),'k-','marker','none');
                hold on;
                for n = 1:size(markers,1)
                    plot(h2,[markers{n,2} markers{n,2}], [0,1], '-','color',markers{n,4},'linewidth',0.5);
                    hold on;
                end
                axis off;
                set(h2,'xlim',[stimOnT(1)-tOffset1 stimOffT(1)+tOffset2]);
            end
            %             set(h1,'ylim',[0 max(maxSpkRealAll(k),maxSpkSponAll)],'xlim',[1 nBins(1,1)]);
            set(h1,'xlim',[1 nBins(1,1)]);
        end
    end



    % extra 2 conditons
    h1 = axes('unit','pixels','pos',[30+200*(5-1) 20+180*(5-1) 180 80]);
    for n = 1:size(markers,1)
        plot(h1,[markers{n,3} markers{n,3}], [0,max(PSTH.maxSpkRealBinMean(k),PSTH.maxSpkSponBinMean)], '--','color',markers{n,4},'linewidth',0.5);
        hold on;
    end
    for r = 1:size((PSTH.spk_data_bin_rate{k,1,5}),2)
        color = ['color',num2str(r)];
        plot(h1,PSTH.spk_data_bin_rate{k,1,5}(:,r)','color',eval(color));
        axis off;
        hold on;
        h2 = axes('unit','pixels','pos',[30+200*(5-1) 20+90+180*(5-1)+8*(r-1) 180 8]);
        stem(h2,stimOnT(1)-tOffset1:stimOffT(1)+tOffset2,spk_data{k,1,5}(stimOnT(1)-tOffset1:stimOffT(1)+tOffset2,r),'k-','marker','none');
        hold on;
        for n = 1:size(markers,1)
            plot(h2,[markers{n,2} markers{n,2}], [0,1], '-','color',markers{n,4},'linewidth',0.5);
            hold on;
        end

        axis off;
        set(h2,'xlim',[stimOnT(1)-tOffset1 stimOffT(1)+tOffset2]);
    end
    set(h1,'xlim',[1 nBins(1,1)]);

    h1 = axes('unit','pixels','pos',[30+200*(5-1) 20+180*(5-5) 180 80]);
    for n = 1:size(markers,1)
        plot(h1,[markers{n,3} markers{n,3}], [0,max(PSTH.maxSpkRealBinMean(k),PSTH.maxSpkSponBinMean)], '--','color',markers{n,4},'linewidth',0.5);
        hold on;
    end
    for r = 1:size((PSTH.spk_data_bin_rate{k,5,5}),2)
        color = ['color',num2str(r)];
        plot(h1,PSTH.spk_data_bin_rate{k,5,5}(:,r)','color',eval(color));
        axis off;
        hold on;
        h2 = axes('unit','pixels','pos',[30+200*(5-1) 20+90+180*(5-5)+8*(r-1) 180 8]);
        stem(h2,stimOnT(1)-tOffset1:stimOffT(1)+tOffset2,spk_data{k,5,5}(stimOnT(1)-tOffset1:stimOffT(1)+tOffset2,r),'k-','marker','none');
        hold on;
        for n = 1:size(markers,1)
            plot(h2,[markers{n,2} markers{n,2}], [0,1], '-','color',markers{n,4},'linewidth',0.5);
            hold on;
        end

        axis off;
        set(h2,'xlim',[stimOnT(1)-tOffset1 stimOffT(1)+tOffset2]);
    end
    set(h1,'xlim',[1 nBins(1,1)]);


    % spontaneous response
    h1 = axes('unit','pixels','pos',[30+200*(1-1) 20+180*(5-1) 180 80]);
    for n = 1:size(markers,1)
        plot(h1,[markers{n,3} markers{n,3}], [0,max(PSTH.maxSpkRealBinMean(k),PSTH.maxSpkSponBinMean)], '--','color',markers{n,4},'linewidth',0.5);
        hold on;
    end
    for r = 1:size(PSTH.spon_spk_data_bin_rate,2)
        color = ['color',num2str(r)];
        plot(h1,PSTH.spon_spk_data_bin_rate(:,r)','color',eval(color));
        axis off;
        hold on;
        h2 = axes('unit','pixels','pos',[30+200*(1-1) 20+90+180*(5-1)+8*(r-1) 180 8]);
        stem(h2,stimOnT(1)-tOffset1:stimOffT(1)+tOffset2,spon_spk_data(stimOnT(1)-tOffset1:stimOffT(1)+tOffset2,r),'k-','marker','none');
        hold on;
        for n = 1:size(markers,1)
            plot(h2,[markers{n,2} markers{n,2}], [0,1], '-','color',markers{n,4},'linewidth',0.5);
            hold on;
        end
        axis off;
        set(h2,'xlim',[stimOnT(1)-tOffset1 stimOffT(1)+tOffset2]);
    end
    set(h1,'xlim',[1 nBins(1,1)]);

    %text on the figure
    axes('unit','pixels','pos',[60 850 1800 80]);
    xlim([0,100]);
    ylim([0,10]);
    if Protocol == DIRECTION_TUNING_3D
        text(30,10,'PSTHs and raster plot for each direction across trials(T)','fontsize',20);
    elseif Protocol == ROTATION_TUNING_3D
        text(30,10,'PSTHs and raster plot for each direction across trials(R)','fontsize',20);
    end
    text(-1,10,'spontaneous','fontsize',15);
    FileNameTemp = num2str(FILE);
    FileNameTemp =  FileNameTemp(1:end);
    str = [FileNameTemp,'_Ch' num2str(SpikeChan)];
    str1 = [FileNameTemp,'\_Ch' num2str(SpikeChan),'    ',stimType{k}];
    text(70,0,str1,'fontsize',18);
    axis off;

    axes('unit','pixels','pos',[60 50 1800 180]);
    xlim([0,100]);
    ylim([0,10]);
    text(0,0,['Max spk FR(Hz): ',num2str(PSTH.maxSpkRealBinAll(k)),' (Bin)'],'fontsize',15);
    text(0,3,['Spon max FR(Hz): ',num2str(maxSpkSponAll), ' (Bin)'],'fontsize',15);

    axis off;

    % save the figure
    str2 = [str,'_PSTHs_raster_trial_',stimType{k}];
    set(gcf,'paperpositionmode','auto');
    switch Protocol
        case DIRECTION_TUNING_3D
            ss = [str2, '_T'];
            saveas(20+k,['Z:\LBY\Recording data\Qiaoqiao\3D_Tuning\Translation\' ss], 'emf');
        case ROTATION_TUNING_3D
            ss = [str2, '_R'];
            saveas(20+k,['Z:\LBY\Recording data\Qiaoqiao\3D_Tuning\Rotation\' ss], 'emf');
    end
end


% % ------ fig.30 plot mean PSTHs across directions (with errorbar)------%
% for k = 1:length(unique_stimType)
%     figure(30+k);
%     set(gcf,'pos',[60 100 1800 900]);
%     clf;
%     [~,h_subplot] = tight_subplot(5,9,0.04,0.15);
% 
%     for j = 2:length(unique_elevation)-1
%         for i = 1:length(unique_azimuth)+1
%             axes(h_subplot(i+(j-1)*9));
%             errorbar(PSTH.spk_data_bin_mean_rate{k}(j,iAzi(i),:),PSTH.spk_data_bin_mean_rate_ste{k}(j,iAzi(i),:),'color','k');
%             hold on;
%             for n = 1:size(markers,1)
%                 plot([markers{n,3} markers{n,3}], [0,max(PSTH.maxSpkRealBinMean(k),PSTH.maxSpkSponBinMean)], '--','color',markers{n,4},'linewidth',0.5);
%                 hold on;
%             end
%             set(gca,'ylim',[0 max(PSTH.maxSpkRealBinMean(k)+PSTH.maxSpkRealBinMeanSte(k),PSTH.maxSpkSponBinMean+PSTH.maxSpkSponBinMeanSte)],'xlim',[1 nBins(1,1)]);
%             %             axis off;
%             SetFigure(15);
%             set(gca,'xtick',[],'xticklabel',[]);
%         end
%     end
% 
%     % 2 extra conditions
%     axes(h_subplot(5+(1-1)*9));
%     errorbar(PSTH.spk_data_bin_mean_rate{k}(1,iAzi(5),:),PSTH.spk_data_bin_mean_rate_ste{k}(1,iAzi(5),:),'color','k');
%     hold on;
%     for n = 1:size(markers,1)
%         plot([markers{n,3} markers{n,3}], [0,max(PSTH.maxSpkRealBinMean(k),PSTH.maxSpkSponBinMean)], '--','color',markers{n,4},'linewidth',0.5);
%         hold on;
%     end
%     set(gca,'ylim',[0 max(PSTH.maxSpkRealBinMean(k)+PSTH.maxSpkRealBinMeanSte(k),PSTH.maxSpkSponBinMean+PSTH.maxSpkSponBinMeanSte)],'xlim',[1 nBins(1,1)]);
%     SetFigure(15);
%     set(gca,'xtick',[],'xticklabel',[]);
% 
%     axes(h_subplot(5+(5-1)*9));
%     errorbar(PSTH.spk_data_bin_mean_rate{k}(5,iAzi(5),:),PSTH.spk_data_bin_mean_rate_ste{k}(5,iAzi(5),:),'color','k');
%     hold on;
%     for n = 1:size(markers,1)
%         plot([markers{n,3} markers{n,3}], [0,max(PSTH.maxSpkRealBinMean(k),PSTH.maxSpkSponBinMean)], '--','color',markers{n,4},'linewidth',0.5);
%         hold on;
%     end
%     set(gca,'ylim',[0 max(PSTH.maxSpkRealBinMean(k)+PSTH.maxSpkRealBinMeanSte(k),PSTH.maxSpkSponBinMean+PSTH.maxSpkSponBinMeanSte)],'xlim',[1 nBins(1,1)]);
% 
%     SetFigure(15);
%     set(gca,'xtick',[],'xticklabel',[]);
% 
%     % spontaneous
%     axes(h_subplot(1+(1-1)*9));
%     errorbar(PSTH.spon_spk_data_bin_mean_rate,PSTH.spon_spk_data_bin_mean_rate_ste,'color','k');
%     hold on;
%     for n = 1:size(markers,1)
%         plot([markers{n,3} markers{n,3}], [0,max(PSTH.maxSpkRealBinMean(k),PSTH.maxSpkSponBinMean)], '--','color',markers{n,4},'linewidth',0.5);
%         hold on;
%     end
%     set(gca,'ylim',[0 max(PSTH.maxSpkRealBinMean(k)+PSTH.maxSpkRealBinMeanSte(k),PSTH.maxSpkSponBinMean+PSTH.maxSpkSponBinMeanSte)],'xlim',[1 nBins(1,1)]);
%     SetFigure(15);
%     set(gca,'xtick',[],'xticklabel',[]);
%     % text on the figure
%     axes('unit','pixels','pos',[60 810 1800 80]);
%     xlim([0,100]);
%     ylim([0,10]);
%     switch Protocol
%         case DIRECTION_TUNING_3D
%             text(36,3,'PSTHs for each direction(T)','fontsize',20);
%         case ROTATION_TUNING_3D
%             text(36,3,'PSTHs for each direction(R)','fontsize',20);
%     end
%     text(1,0,'spontaneous','fontsize',15);
%     FileNameTemp = num2str(FILE);
%     FileNameTemp =  FileNameTemp(1:end);
%     str = [FileNameTemp,'_Ch' num2str(SpikeChan)];
%     str1 = [FileNameTemp,'\_Ch' num2str(SpikeChan),'    ',stimType{k}];
%     text(70,0,str1,'fontsize',18);
%     axis off;
% 
%     axes('unit','pixels','pos',[60 80 1800 100]);
%     xlim([0,100]);
%     ylim([0,10]);
%     text(0,7,['Spon max FR(Hz): ',num2str(PSTH.maxSpkSponBinMean), ' (Bin)'],'fontsize',15);
%     text(0,10,['Spon mean FR(Hz): ',num2str(PSTH.meanSpkSponBinMean), ' (Bin)'],'fontsize',15);
% 
% 
%     %--- this is for annotation - the direction ----%
%     text(4,0,'\downarrow ','fontsize',30);
%     text(24,0,'\leftarrow ','fontsize',30);
%     text(46,0,'\uparrow ','fontsize',30);
%     text(67,0,'\rightarrow ','fontsize',30);
%     text(88,0,'\downarrow ','fontsize',30);
%     text(52,65,'\uparrow ','fontsize',30);
%     text(52,10,'\downarrow ','fontsize',30);
%     %--- this is for annotation - the direction ----%
%     axis off;
% 
%     % to save the figures
%     str3 = [str, '_PSTH_',stimType{k}];
%     set(gcf,'paperpositionmode','auto');
%     if Protocol == DIRECTION_TUNING_3D
%         ss = [str3, '_T'];
%         saveas(30+k,['Z:\LBY\Recording data\Qiaoqiao\3D_Tuning\Translation\' ss], 'emf');
%     elseif Protocol == ROTATION_TUNING_3D
%         ss = [str3, '_R'];
%         saveas(30+k,['Z:\LBY\Recording data\Qiaoqiao\3D_Tuning\Rotation\' ss], 'emf');
%     end
% 
% end

% ------ fig.40 plot mean PSTHs across directions (normal/with peak time)------%

peakColor = {colorDRed,colorDBlue,colorDGreen,colorLRed,colorLBlue,colorDGray,colorDOrange};

for k = 1:length(unique_stimType)
% for k = 1
    figure(40+k);
    set(gcf,'pos',[60 100 1800 900]);
    clf;
    [~,h_subplot] = tight_subplot(5,9,0.04,0.15);
    
    for j = 2:length(unique_elevation)-1
        for i = 1:length(unique_azimuth)+1
            axes(h_subplot(i+(j-1)*9));
            plot(squeeze(PSTH.spk_data_bin_mean_rate{k}(j,iAzi(i),:)),'color','k','linewidth',2);
            hold on;
%                         for n = 1:length(PSTH.peak{k})
%                             plot([PSTH.peak{k}(n) PSTH.peak{k}(n)], [0,max(PSTH.maxSpkRealBinMean(k),PSTH.maxSpkSponBinMean)], '--','color',peakColor{n},'linewidth',1.5);
%                             hold on;
%                         end
            for n = 1:size(markers,1)
                plot([markers{n,3} markers{n,3}], [0,max(PSTH.maxSpkRealBinMean(k),PSTH.maxSpkSponBinMean)], '--','color',markers{n,4},'linewidth',0.5);
                hold on;
            end
            set(gca,'ylim',[0 max(PSTH.maxSpkRealBinMean(k)+PSTH.maxSpkRealBinMeanSte(k),PSTH.maxSpkSponBinMean+PSTH.maxSpkSponBinMeanSte)],'xlim',[1 nBins(1,1)]);
            %             axis off;
            SetFigure(15);
            set(gca,'xtick',[],'xticklabel',[]);
        end
    end
    
    % 2 extra conditions
    axes(h_subplot(5+(1-1)*9));
    plot(squeeze(PSTH.spk_data_bin_mean_rate{k}(1,iAzi(5),:)),'color','k','linewidth',2);
    hold on;
%         for n = 1:length(PSTH.peak{k})
%             plot([PSTH.peak{k}(n) PSTH.peak{k}(n)], [0,max(PSTH.maxSpkRealBinMean(k),PSTH.maxSpkSponBinMean)], '--','color',peakColor{n},'linewidth',1.5);
%             hold on;
%         end
    for n = 1:size(markers,1)
        plot([markers{n,3} markers{n,3}], [0,max(PSTH.maxSpkRealBinMean(k),PSTH.maxSpkSponBinMean)], '--','color',markers{n,4},'linewidth',0.5);
        hold on;
    end
    set(gca,'ylim',[0 max(PSTH.maxSpkRealBinMean(k)+PSTH.maxSpkRealBinMeanSte(k),PSTH.maxSpkSponBinMean+PSTH.maxSpkSponBinMeanSte)],'xlim',[1 nBins(1,1)]);
    SetFigure(15);
    set(gca,'xtick',[],'xticklabel',[]);
    
    axes(h_subplot(5+(5-1)*9));
    plot(squeeze(PSTH.spk_data_bin_mean_rate{k}(5,iAzi(5),:)),'color','k','linewidth',2);
    hold on;
%         for n = 1:length(PSTH.peak{k})
%             plot([PSTH.peak{k}(n) PSTH.peak{k}(n)], [0,max(PSTH.maxSpkRealBinMean(k),PSTH.maxSpkSponBinMean)], '--','color',peakColor{n},'linewidth',1.5);
%             hold on;
%         end
    for n = 1:size(markers,1)
        plot([markers{n,3} markers{n,3}], [0,max(PSTH.maxSpkRealBinMean(k),PSTH.maxSpkSponBinMean)], '--','color',markers{n,4},'linewidth',0.5);
        hold on;
    end
    set(gca,'ylim',[0 max(PSTH.maxSpkRealBinMean(k)+PSTH.maxSpkRealBinMeanSte(k),PSTH.maxSpkSponBinMean+PSTH.maxSpkSponBinMeanSte)],'xlim',[1 nBins(1,1)]);
    
    SetFigure(15);
    set(gca,'xtick',[],'xticklabel',[]);
    
    % spontaneous
    axes(h_subplot(1+(1-1)*9));
    plot(squeeze(PSTH.spon_spk_data_bin_mean_rate),'color','k','linewidth',2);
    hold on;
    
    set(gca,'ylim',[0 max(PSTH.maxSpkRealBinMean(k)+PSTH.maxSpkRealBinMeanSte(k),PSTH.maxSpkSponBinMean+PSTH.maxSpkSponBinMeanSte)],'xlim',[1 nBins(1,1)]);
    SetFigure(15);
    set(gca,'xtick',[],'xticklabel',[]);
    % text on the figure
    axes('unit','pixels','pos',[60 810 1800 80]);
    xlim([0,100]);
    ylim([0,10]);
    switch Protocol
        case DIRECTION_TUNING_3D
            text(36,3,'PSTHs for each direction(T)','fontsize',20);
        case ROTATION_TUNING_3D
            text(36,3,'PSTHs for each direction(R)','fontsize',20);
    end
    text(1,0,'Spontaneous','fontsize',15);
    FileNameTemp = num2str(FILE);
    FileNameTemp =  FileNameTemp(1:end);
    str = [FileNameTemp,'_Ch' num2str(SpikeChan)];
    str1 = [FileNameTemp,'\_Ch' num2str(SpikeChan),'    ',stimType{k}];
    text(70,0,str1,'fontsize',18);
%     for n = 1:length(PSTH.peak{k})
%         text(60,-5*n,['t',num2str(n),' = ',num2str((PSTH.peak{k}(n)*timeStep-tOffset1)/1000),' s'],'fontsize',15);
%     end
    axis off;
    
    axes('unit','pixels','pos',[60 80 1800 100]);
    xlim([0,100]);
    ylim([0,10]);
    text(0,7,['Spon max FR(Hz): ',num2str(PSTH.maxSpkSponBinMean), ' (Bin)'],'fontsize',15);
    text(0,10,['Spon mean FR(Hz): ',num2str(PSTH.meanSpkSponBinMean), ' (Bin)'],'fontsize',15);
    
    
%     %--- this is for annotation - the direction ----%
%     text(4,0,'\downarrow ','fontsize',30);
%     text(24,0,'\leftarrow ','fontsize',30);
%     text(46,0,'\uparrow ','fontsize',30);
%     text(67,0,'\rightarrow ','fontsize',30);
%     text(88,0,'\downarrow ','fontsize',30);
%     text(52,65,'\uparrow ','fontsize',30);
%     text(52,10,'\downarrow ','fontsize',30);
%     %--- this is for annotation - the direction ----%
    axis off;
    
    % to save the figures
    str3 = [str, '_PSTH_peak',stimType{k}];
    set(gcf,'paperpositionmode','auto');
    if Protocol == DIRECTION_TUNING_3D
        ss = [str3, '_T'];
        saveas(40+k,['Z:\LBY\Recording data\Qiaoqiao\3D_Tuning\Translation\' ss], 'emf');
    elseif Protocol == ROTATION_TUNING_3D
        ss = [str3, '_R'];
        saveas(40+k,['Z:\LBY\Recording data\Qiaoqiao\3D_Tuning\Rotation\' ss], 'emf');
    end
    
end

% % ------ fig.60 plot mean PSTHs across directions(with DS peak time) ------%
% 
% peakColor = {colorDRed,colorDBlue,colorDGreen,colorLRed};
% 
% for k = 1:length(unique_stimType)
%     figure(60+k);
%     set(gcf,'pos',[60 100 1800 900]);
%     clf;
%     [~,h_subplot] = tight_subplot(5,9,0.04,0.15);
%     
%     for j = 2:length(unique_elevation)-1
%         for i = 1:length(unique_azimuth)+1
%             axes(h_subplot(i+(j-1)*9));
%             plot(squeeze(PSTH.spk_data_bin_mean_rate{k}(j,iAzi(i),:)),'color','k','linewidth',2);
%             hold on;
%             for n = 1:length(PSTH.peak_DS{k})
%                 plot([PSTH.peak_DS{k}(n) PSTH.peak_DS{k}(n)], [0,max(PSTH.maxSpkRealBinMean(k),PSTH.maxSpkSponBinMean)], '--','color',peakColor{n},'linewidth',1.5);
%                 hold on;
%             end
%             set(gca,'ylim',[0 max(PSTH.maxSpkRealBinMean(k)+PSTH.maxSpkRealBinMeanSte(k),PSTH.maxSpkSponBinMean+PSTH.maxSpkSponBinMeanSte)],'xlim',[1 nBins(1,1)]);
%             %             axis off;
%             SetFigure(15);
%             set(gca,'xtick',[],'xticklabel',[]);
%         end
%     end
%     
%     % 2 extra conditions
%     axes(h_subplot(5+(1-1)*9));
%     plot(squeeze(PSTH.spk_data_bin_mean_rate{k}(1,iAzi(5),:)),'color','k','linewidth',2);
%     hold on;
%     for n = 1:length(PSTH.peak_DS{k})
%         plot([PSTH.peak_DS{k}(n) PSTH.peak_DS{k}(n)], [0,max(PSTH.maxSpkRealBinMean(k),PSTH.maxSpkSponBinMean)], '--','color',peakColor{n},'linewidth',1.5);
%         hold on;
%     end
%     set(gca,'ylim',[0 max(PSTH.maxSpkRealBinMean(k)+PSTH.maxSpkRealBinMeanSte(k),PSTH.maxSpkSponBinMean+PSTH.maxSpkSponBinMeanSte)],'xlim',[1 nBins(1,1)]);
%     SetFigure(15);
%     set(gca,'xtick',[],'xticklabel',[]);
%     
%     axes(h_subplot(5+(5-1)*9));
%     plot(squeeze(PSTH.spk_data_bin_mean_rate{k}(5,iAzi(5),:)),'color','k','linewidth',2);
%     hold on;
%     for n = 1:length(PSTH.peak_DS{k})
%         plot([PSTH.peak_DS{k}(n) PSTH.peak_DS{k}(n)], [0,max(PSTH.maxSpkRealBinMean(k),PSTH.maxSpkSponBinMean)], '--','color',peakColor{n},'linewidth',1.5);
%         hold on;
%     end
%     set(gca,'ylim',[0 max(PSTH.maxSpkRealBinMean(k)+PSTH.maxSpkRealBinMeanSte(k),PSTH.maxSpkSponBinMean+PSTH.maxSpkSponBinMeanSte)],'xlim',[1 nBins(1,1)]);
%     
%     SetFigure(15);
%     set(gca,'xtick',[],'xticklabel',[]);
%     
%     % spontaneous
%     axes(h_subplot(1+(1-1)*9));
%     plot(squeeze(PSTH.spon_spk_data_bin_mean_rate),'color','k','linewidth',2);
%     hold on;
%     
%     set(gca,'ylim',[0 max(PSTH.maxSpkRealBinMean(k)+PSTH.maxSpkRealBinMeanSte(k),PSTH.maxSpkSponBinMean+PSTH.maxSpkSponBinMeanSte)],'xlim',[1 nBins(1,1)]);
%     SetFigure(15);
%     set(gca,'xtick',[],'xticklabel',[]);
%     % text on the figure
%     axes('unit','pixels','pos',[60 810 1800 80]);
%     xlim([0,100]);
%     ylim([0,10]);
%     switch Protocol
%         case DIRECTION_TUNING_3D
%             text(36,3,'PSTHs for each direction(T)','fontsize',20);
%         case ROTATION_TUNING_3D
%             text(36,3,'PSTHs for each direction(R)','fontsize',20);
%     end
%     text(1,0,'spontaneous','fontsize',15);
%     FileNameTemp = num2str(FILE);
%     FileNameTemp =  FileNameTemp(1:end);
%     str = [FileNameTemp,'_Ch' num2str(SpikeChan)];
%     str1 = [FileNameTemp,'\_Ch' num2str(SpikeChan),'    ',stimType{k}];
%     text(70,0,str1,'fontsize',18);
%     for n = 1:length(PSTH.peak_DS{k})
%         text(60,-5*n,['t',num2str(n),' = ',num2str((PSTH.peak_DS{k}(n)*timeStep-tOffset1)/1000),' s'],'fontsize',15);
%     end
%     axis off;
%     
%     axes('unit','pixels','pos',[60 80 1800 100]);
%     xlim([0,100]);
%     ylim([0,10]);
%     text(0,7,['Spon max FR(Hz): ',num2str(PSTH.maxSpkSponBinMean), ' (Bin)'],'fontsize',15);
%     text(0,10,['Spon mean FR(Hz): ',num2str(PSTH.meanSpkSponBinMean), ' (Bin)'],'fontsize',15);
%     
%     
% %     %--- this is for annotation - the direction ----%
% %     text(4,0,'\downarrow ','fontsize',30);
% %     text(24,0,'\leftarrow ','fontsize',30);
% %     text(46,0,'\uparrow ','fontsize',30);
% %     text(67,0,'\rightarrow ','fontsize',30);
% %     text(88,0,'\downarrow ','fontsize',30);
% %     text(52,65,'\uparrow ','fontsize',30);
% %     text(52,10,'\downarrow ','fontsize',30);
% %     %--- this is for annotation - the direction ----%
%     axis off;
%     
%     % to save the figures
%     str3 = [str, '_PSTH_peak_DS',stimType{k}];
%     set(gcf,'paperpositionmode','auto');
%     if Protocol == DIRECTION_TUNING_3D
%         ss = [str3, '_T'];
%         saveas(60+k,['Z:\LBY\Recording data\Qiaoqiao\3D_Tuning\Translation\' ss], 'emf');
%     elseif Protocol == ROTATION_TUNING_3D
%         ss = [str3, '_R'];
%         saveas(60+k,['Z:\LBY\Recording data\Qiaoqiao\3D_Tuning\Rotation\' ss], 'emf');
%     end
%     
% end


% % ------ fig.50 plot mean PSTHs-spon across directions ------%
% 
% 
% for k = 1:length(unique_stimType)
%     figure(50+k);
%     set(gcf,'pos',[60 100 1800 900]);
%     clf;
%     maxY = max(PSTH.maxSpkRealBinMean(k)-PSTH.spon_spk_data_bin_mean_rate);
%     minY = min(PSTH.minSpkRealBinMean(k)-PSTH.spon_spk_data_bin_mean_rate);
%     for j = 2:length(unique_elevation)-1
%         for i = 1:length(unique_azimuth)+1
%             subplot(5,9,i+(j-1)*9);
%             plot(squeeze(PSTH.spk_data_bin_mean_rate{k}(j,iAzi(i),:))-PSTH.spon_spk_data_bin_mean_rate,'color','k','linewidth',1.5);
%             hold on;
%             plot([1 nBins(1,1)],[0 0],':k','linewidth',0.5);
%             hold on;
%             for n = 1:size(markers,1)
%                 plot([markers{n,3} markers{n,3}],[minY maxY], '--','color',markers{n,4},'linewidth',0.5);
%                 hold on;
%             end
%             set(gca,'ylim',[minY maxY],'xlim',[1 nBins(1,1)]);
%             axis off;
%         end
%     end
%     % 2 extra conditions
%     subplot(5,9,5+(1-1)*9);
%     plot(squeeze(PSTH.spk_data_bin_mean_rate{k}(1,iAzi(5),:))-PSTH.spon_spk_data_bin_mean_rate,'color','k','linewidth',1.5);
%     hold on;
%     plot([1 nBins(1,1)],[0 0],':k','linewidth',0.5);
%     hold on;
%     for n = 1:size(markers,1)
%         plot([markers{n,3} markers{n,3}],[minY maxY], '--','color',markers{n,4},'linewidth',0.5);
%         hold on;
%     end
%     set(gca,'ylim',[minY maxY],'xlim',[1 nBins(1,1)]);
%     axis off;
% 
%     subplot(5,9,5+(5-1)*9);
%     plot(squeeze(PSTH.spk_data_bin_mean_rate{k}(5,iAzi(5),:))-PSTH.spon_spk_data_bin_mean_rate,'color','k','linewidth',1.5);
%     hold on;
%     plot([1 nBins(1,1)],[0 0],':k','linewidth',0.5);
%     hold on;
%     for n = 1:size(markers,1)
%         plot([markers{n,3} markers{n,3}], [minY maxY], '--','color',markers{n,4},'linewidth',0.5);
%         hold on;
%     end
%     set(gca,'ylim',[minY maxY],'xlim',[1 nBins(1,1)]);
%     axis off;
% 
%     % spontaneous
%     subplot(5,9,1+(1-1)*9);
%     errorbar(PSTH.spon_spk_data_bin_mean_rate,PSTH.spon_spk_data_bin_mean_rate_ste,'color','k');
%     hold on;
%     for n = 1:size(markers,1)
%         plot([markers{n,3} markers{n,3}], [0,max(PSTH.maxSpkRealBinMean(k),PSTH.maxSpkSponBinMean)], '--','color',markers{n,4},'linewidth',0.5);
%         hold on;
%     end
%     set(gca,'ylim',[0 max(PSTH.maxSpkRealBinMean(k)+PSTH.maxSpkRealBinMeanSte(k),PSTH.maxSpkSponBinMean+PSTH.maxSpkSponBinMeanSte)],'xlim',[1 nBins(1,1)]);
%     axis off;
% 
% 
%     % text on the figure
%     axes('unit','pixels','pos',[60 830 1800 80]);
%     xlim([0,100]);
%     ylim([0,10]);
%     text(40,3,'PSTHs for each direction (-spon)','fontsize',20);
%     text(11,2,'spontaneous','fontsize',15);
%     FileNameTemp = num2str(FILE);
%     FileNameTemp =  FileNameTemp(1:end);
%     str = [FileNameTemp,'_Ch' num2str(SpikeChan)];
%     str1 = [FileNameTemp,'\_Ch' num2str(SpikeChan),'    ',stimType{k}];
%     text(70,0,str1,'fontsize',18);
%     axis off;
% 
%     axes('unit','pixels','pos',[60 50 1800 100]);
%     xlim([0,100]);
%     ylim([0,10]);
%     %     text(0,0,['Max spk mean FR(Hz): ',num2str(maxSpkRealAll(k))],'fontsize',15);
%     text(0,2,['Spon max mean FR(Hz): ',num2str(PSTH.maxSpkSponBinMean)],'fontsize',15);
%     text(0,4,['Spon mean FR(Hz): ',num2str(PSTH.meanSpkSponBinMean)],'fontsize',15);
%     axis off;
% 
% %     % to save the figures
% %     str3 = [str, '_PSTH-spon_',stimType{k}];
% %     set(gcf,'paperpositionmode','auto');
% %     if Protocol == DIRECTION_TUNING_3D
% %         ss = [str3, '_T'];
% %         saveas(50+k,['Z:\LBY\Recording data\Qiaoqiao\3D_Tuning\Translation\' ss], 'emf');
% %     elseif Protocol == ROTATION_TUNING_3D
% %         ss = [str3, '_R'];
% %         saveas(50+k,['Z:\LBY\Recording data\Qiaoqiao\3D_Tuning\Rotation\' ss], 'emf');
% %     end

% end
