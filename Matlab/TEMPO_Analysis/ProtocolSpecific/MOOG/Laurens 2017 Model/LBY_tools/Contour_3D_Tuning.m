%% plot figures - contour
%

% % ------ fig.110 plot coutour tuning responses (sum response of total 1.5s/ middle ?) ------%
% 
% % transform
% % k=1,2,3
% % j= -90,-45,0,45,90 (up->down)
% % i=0 45 90 135 180 225 270 315
% % 270-225-180-135-90-45-0-315-270 for figures
% iAzi = [7 6 5 4 3 2 1 8 7];
% for k = 1:length(unique_stimType)
%     spk_data_count_mean_rate_trans{k} = [];
%     for i = 1:length(unique_azimuth)+1
%         spk_data_count_mean_rate_trans{k}(:,i) = PSTH.spk_data_count_mean_rate{k}(:,iAzi(i));
%     end
% end
% xAzi = 1:9;
% yEle = [-1,-0.707,0,0.707,1];
% 
% figure(110);
% set(gcf,'pos',[60 200 1800 800]);
% clf;
% 
% for k = 1:length(unique_stimType)
%     axes('unit','pixels','pos',[60+100*k+700*(k-1) 300 650 300]);
%     contourf(xAzi,yEle,spk_data_count_mean_rate_trans{k},'linecolor','w','linestyle','none','linewidth',3);
%     a = cell2mat(spk_data_count_mean_rate_trans);
%     caxis([min(a(:)) max(a(:))]);
%     colorbar;
%     set(gca, 'ydir' , 'reverse'); % so that up is up, down is down
%     title(gca,stimType{k});
%     set(gca, 'xtick', [] );
%     set(gca, 'ytick', [] );
%     box off;
%     axes('unit','pixels','pos',[60+100*k+700*(k-1) 200 567 100]);
%     for i=1:(length(unique_azimuth)+1 )
%         y_azimuth_mean(1,i)=mean(spk_data_count_mean_rate_trans{k}(:,i));
%         y_azimuth_std(1,i) =std( spk_data_count_mean_rate_trans{k}(:,i));
%         y_azimuth_ste(1,i) =y_azimuth_std(1,i)/ sqrt(length(find( (temp_azimuth==unique_azimuth(iAzi(i)))&(temp_stimType==unique_stimType(k)) )) );
%     end
%     % errobar of azimuth
%     errorbar(xAzi,y_azimuth_mean,y_azimuth_ste,'k.-','markerSize',20,'linewidth',2);
%     xlim([1 9]);
%     set(gca, 'XTickMode','manual');
%     set(gca, 'xtick',[1,2,3,4,5,6,7,8,9]);
%     set(gca, 'xticklabel','270|225|180|135|90|45|0|-45|-90');
%     xlabel('Azimuth');
%     ylim([min(y_azimuth_mean(1,:))-max(y_azimuth_ste(1,:)), max(y_azimuth_mean(1,:))+max(y_azimuth_ste(1,:))+2]);
%     axis off;
%     % errobar of elevation
%     axes('unit','pixels','pos',[60+800*(k-1) 300 100 300]);
%     for j=1:length(unique_elevation)
%         y_elevation_mean(1,j)=mean(spk_data_count_mean_rate_trans{k}(j,:));
%         y_elevation_std(1,j) =std( spk_data_count_mean_rate_trans{k}(j,:));
%         y_elevation_ste(1,j) =y_elevation_std(1,j)/ sqrt(length(find( (temp_elevation==unique_elevation(j))&(temp_stimType==unique_stimType(k)) )) );
%     end
%     errorbar(yEle,y_elevation_mean,y_elevation_ste,'k.-','markerSize',20,'linewidth',2);
%     xlim([-1 1]);
%     view(90,270);
%     xlabel('Elevation');
%     ylim([min(y_elevation_mean(1,:))-max(y_elevation_ste(1,:)), max(y_elevation_mean(1,:))+max(y_elevation_ste(1,:))+2]);
%     axis off;
% 
%     % text on the figure
%     axes('unit','pixels','pos',[60 730 1800 80]);
%     xlim([0,100]);
%     ylim([0,10]);
%     text(25,0,'Color contour maps','fontsize',28);
%     FileNameTemp = num2str(FILE);
%     FileNameTemp =  FileNameTemp(1:end);
%     str = [FileNameTemp,'_Ch' num2str(SpikeChan)];
%     str1 = [FileNameTemp,'\_Ch' num2str(SpikeChan)];
%     text(60,0,str1,'fontsize',24);
%     axis off;
% 
%     axes('unit','pixels','pos',[60 50 1800 100]);
%     xlim([0,100]);
%     ylim([0,10]);
%     text(0,0,'Max spk FR(Hz): ','fontsize',15);
%     text(11+11*(k-1),0,[stimType{k},': ',num2str(maxSpkRealMean(k))],'fontsize',13);
%     text(0,3,['Spon max FR(Hz): ',num2str(maxSpkSponMean)],'fontsize',15);
%     text(0,6,['Spon mean FR(Hz): ',num2str(meanSpkSponMean)],'fontsize',15);
%     axis off;
% 
%     str4 = [str, '_Tuning_totalT'];
%     if Protocol == DIRECTION_TUNING_3D
%         ss = [str4, '_T'];
%         saveas(gcf,['Z:\LBY\Recording data\Qiaoqiao\3D_Tuning\Translation\' ss], 'emf');
%     elseif Protocol == ROTATION_TUNING_3D
%         ss = [str4, '_R'];
%         saveas(gcf,['Z:\LBY\Recording data\Qiaoqiao\3D_Tuning\Rotation\' ss], 'emf');
%     end
% 
% end
% 
% % ------ fig.120 plot coutour tuning responses (at peak time) ------%
% 
% % transform
% % k=1,2,3
% % j= -90,-45,0,45,90 (up->down)
% % i=0 45 90 135 180 225 270 315
% % 270-225-180-135-90-45-0-315-270 for figures
% iAzi = [7 6 5 4 3 2 1 8 7];
% for k = 1:length(unique_stimType)
%     spk_data_bin_mean_rate_trans{k} = [];
%     for i = 1:length(unique_azimuth)+1
%         spk_data_bin_mean_rate_trans{k}(:,i,:) = PSTH.spk_data_bin_mean_rate{k}(:,iAzi(i),:);
%     end
% end
% xAzi = 1:9;
% yEle = [-1,-0.707,0,0.707,1];
% 
% 
% 
% for k = 1:length(unique_stimType)
%     if ~isempty(PSTH.peak_DS{k})
%         figure(120+k);
%         set(gcf,'pos',[60 200 1800 800]);
%         clf;
%         for ii = 1:length(PSTH.peak_DS{k})
%             axes('unit','pixels','pos',[60+100*ii+500*(ii-1) 300 400 200]);
%             contourf(xAzi,yEle,spk_data_bin_mean_rate_trans{k}(:,:,PSTH.peak_DS{k}(ii)),'linecolor','w','linestyle','none');
%             colorbar;
%             set(gca, 'ydir' , 'reverse'); % so that up is up, down is down
%             title(gca,['t = ',num2str((PSTH.peak_DS{k}(ii)*timeStep-tOffset1)/1000),' s']);
%             set(gca, 'xtick', [] );
%             set(gca, 'ytick', [] );
%             box off;
%             % errobar of azimuth
%             axes('unit','pixels','pos',[60+100*ii+500*(ii-1) 200 315 100]);
%             y_azimuth_mean = mean(spk_data_bin_mean_rate_trans{k}(:,:,PSTH.peak_DS{k}(ii)),1);
%             y_azimuth_std =std(spk_data_bin_mean_rate_trans{k}(:,:,PSTH.peak_DS{k}(ii)),1);
%             y_azimuth_ste =y_azimuth_std / sqrt(length(find( (temp_azimuth==unique_azimuth(iAzi(i)))&(temp_stimType==unique_stimType(k)) )) );
%             
%             errorbar(xAzi,y_azimuth_mean,y_azimuth_ste,'k.-','markerSize',20,'linewidth',2);
%             xlim([1 9]);
%             set(gca, 'XTickMode','manual');
%             set(gca, 'xtick',[1,2,3,4,5,6,7,8,9]);
%             set(gca, 'xticklabel','270|225|180|135|90|45|0|-45|-90');
%             xlabel('Azimuth');
%             ylim([min(y_azimuth_mean)-max(y_azimuth_ste), max(y_azimuth_mean)+max(y_azimuth_ste)+2]);
%             axis off;
%             % errobar of elevation
%             axes('unit','pixels','pos',[60+600*(ii-1) 300 100 200]);
%             y_elevation_mean=mean(spk_data_bin_mean_rate_trans{k}(:,:,PSTH.peak_DS{k}(ii)),2);
%             y_elevation_std =std(spk_data_bin_mean_rate_trans{k}(:,:,PSTH.peak_DS{k}(ii)),0,2);
%             y_elevation_ste =y_elevation_std/ sqrt(length(find( (temp_elevation==unique_elevation(j))&(temp_stimType==unique_stimType(k)) )) );
%             errorbar(yEle,y_elevation_mean,y_elevation_ste,'k.-','markerSize',20,'linewidth',2);
%             xlim([-1 1]);
%             view(90,270);
%             xlabel('Elevation');
%             ylim([min(y_elevation_mean)-max(y_elevation_ste), max(y_elevation_mean)+max(y_elevation_ste)+2]);
%             axis off;
%         end
%         % text on the figure
%         axes('unit','pixels','pos',[60 730 1800 80]);
%         xlim([0,100]);
%         ylim([0,10]);
%         text(25,0,'Color contour maps','fontsize',28);
%         FileNameTemp = num2str(FILE);
%         FileNameTemp =  FileNameTemp(1:end);
%         str = [FileNameTemp,'_Ch' num2str(SpikeChan),'_',stimType{k}];
%         str1 = [FileNameTemp,'\_Ch' num2str(SpikeChan),'    ',stimType{k}];
%         text(60,0,str1,'fontsize',24);
%         axis off;
%         
%         axes('unit','pixels','pos',[60 50 1800 100]);
%         xlim([0,100]);
%         ylim([0,10]);
%         text(0,0,'Max spk FR(Hz): ','fontsize',15);
%         text(11,0,num2str(maxSpkRealMean(k)),'fontsize',13);
%         text(0,3,['Spon max FR(Hz): ',num2str(maxSpkSponMean)],'fontsize',15);
%         text(0,6,['Spon mean FR(Hz): ',num2str(meanSpkSponMean)],'fontsize',15);
%         axis off;
%         
%         if Protocol == DIRECTION_TUNING_3D
%             ss = [str, '_T'];
%             saveas(gcf,['Z:\LBY\Recording data\Qiaoqiao\3D_Tuning\Translation\' ss], 'emf');
%         elseif Protocol == ROTATION_TUNING_3D
%             ss = [str, '_R'];
%             saveas(gcf,['Z:\LBY\Recording data\Qiaoqiao\3D_Tuning\Rotation\' ss], 'emf');
%         end
%     end
% end

% ------ fig.130 plot coutour tuning responses (3 time) ------%

% transform
% k=1,2,3
% j= -90,-45,0,45,90 (up->down)
% i=0 45 90 135 180 225 270 315
% 270-225-180-135-90-45-0-315-270 for figures
iAzi = [7 6 5 4 3 2 1 8 7];
for k = 1:length(unique_stimType)
    spk_data_bin_mean_rate_trans{k} = [];
    for i = 1:length(unique_azimuth)+1
        spk_data_bin_mean_rate_trans{k}(:,i,:) = PSTH.spk_data_bin_mean_rate{k}(:,iAzi(i),:);
    end
end
xAzi = 1:9;
yEle = [-1,-0.707,0,0.707,1];

aMaxBin = round((stimOnT(1)+aMax-PSTH_onT+timeStep)/timeStep);
aMinBin = round((stimOnT(1)+aMin-PSTH_onT+timeStep)/timeStep);
vMaxBin = round((stimOnT(1)+aMax+(aMin-aMax)/2-PSTH_onT+timeStep)/timeStep);

Tmarker = [aMaxBin,vMaxBin,aMinBin];

for k = 1:length(unique_stimType)
    if ~isempty(Tmarker)
        figure(130+k);
        set(gcf,'pos',[60 200 1800 800]);
        clf;
        for ii = 1:3
            axes('unit','pixels','pos',[60+100*ii+500*(ii-1) 300 400 200]);
            contourf(xAzi,yEle,spk_data_bin_mean_rate_trans{k}(:,:,Tmarker(ii)),'linecolor','w','linestyle','none');
                a = spk_data_bin_mean_rate_trans{k};
    caxis([min(a(:)) max(a(:))]);
            colorbar;
            set(gca, 'ydir' , 'reverse'); % so that up is up, down is down
            title(gca,['t = ',num2str((Tmarker(ii)*timeStep-tOffset1)/1000),' s']);
            set(gca, 'xtick', [] );
            set(gca, 'ytick', [] );
            box off;
            % errobar of azimuth
            axes('unit','pixels','pos',[60+100*ii+500*(ii-1) 200 315 100]);
            y_azimuth_mean = mean(spk_data_bin_mean_rate_trans{k}(:,:,Tmarker(ii)),1);
            y_azimuth_std =std(spk_data_bin_mean_rate_trans{k}(:,:,Tmarker(ii)),1);
            y_azimuth_ste =y_azimuth_std / sqrt(length(find( (temp_azimuth==unique_azimuth(iAzi(i)))&(temp_stimType==unique_stimType(k)) )) );
            
            errorbar(xAzi,y_azimuth_mean,y_azimuth_ste,'k.-','markerSize',20,'linewidth',2);
            xlim([1 9]);
            set(gca, 'XTickMode','manual');
            set(gca, 'xtick',[1,2,3,4,5,6,7,8,9]);
            set(gca, 'xticklabel','270|225|180|135|90|45|0|-45|-90');
            xlabel('Azimuth');
            ylim([min(y_azimuth_mean)-max(y_azimuth_ste), max(y_azimuth_mean)+max(y_azimuth_ste)+2]);
            axis off;
            % errobar of elevation
            axes('unit','pixels','pos',[60+600*(ii-1) 300 100 200]);
            y_elevation_mean=mean(spk_data_bin_mean_rate_trans{k}(:,:,Tmarker(ii)),2);
            y_elevation_std =std(spk_data_bin_mean_rate_trans{k}(:,:,Tmarker(ii)),0,2);
            y_elevation_ste =y_elevation_std/ sqrt(length(find( (temp_elevation==unique_elevation(j))&(temp_stimType==unique_stimType(k)) )) );
            errorbar(yEle,y_elevation_mean,y_elevation_ste,'k.-','markerSize',20,'linewidth',2);
            xlim([-1 1]);
            view(90,270);
            xlabel('Elevation');
            ylim([min(y_elevation_mean)-max(y_elevation_ste), max(y_elevation_mean)+max(y_elevation_ste)+2]);
            axis off;
        end
        % text on the figure
        axes('unit','pixels','pos',[60 730 1800 80]);
        xlim([0,100]);
        ylim([0,10]);
        text(25,0,'Color contour maps','fontsize',28);
        FileNameTemp = num2str(FILE);
        FileNameTemp =  FileNameTemp(1:end);
        str = [FileNameTemp,'_Ch' num2str(SpikeChan),'_',stimType{k}];
        str1 = [FileNameTemp,'\_Ch' num2str(SpikeChan),'    ',stimType{k}];
        text(60,0,str1,'fontsize',24);
        axis off;
        
        axes('unit','pixels','pos',[60 50 1800 100]);
        xlim([0,100]);
        ylim([0,10]);
        text(0,0,'Max spk FR(Hz): ','fontsize',15);
        text(11,0,num2str(maxSpkRealMean(k)),'fontsize',13);
        text(0,3,['Spon max FR(Hz): ',num2str(maxSpkSponMean)],'fontsize',15);
        text(0,6,['Spon mean FR(Hz): ',num2str(meanSpkSponMean)],'fontsize',15);
        axis off;
        
        if Protocol == DIRECTION_TUNING_3D
            ss = [str, '_T'];
            saveas(gcf,['Z:\LBY\Recording data\Qiaoqiao\3D_Tuning\Translation\' ss], 'emf');
        elseif Protocol == ROTATION_TUNING_3D
            ss = [str, '_R'];
            saveas(gcf,['Z:\LBY\Recording data\Qiaoqiao\3D_Tuning\Rotation\' ss], 'emf');
        end
    end
end

