% plot ontour figures across time



% ------ fig.130 plot coutour tuning responses (across time) ------%

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


% from stimOnBin to stimOffBin

for k = 1:length(unique_stimType)
    pic_num = 1;
    for ii = stimOnBin : stimOffBin
        clf;
        figure(130+k);
        set(gcf,'pos',[60 200 1000 800],'color','w');
        axes('pos',[0.2 0.2 0.7 0.5]);
        contourf(xAzi,yEle,spk_data_bin_mean_rate_trans{k}(:,:,ii),'linecolor','w','linestyle','none');
        a = spk_data_bin_mean_rate_trans{k};
        caxis([min(a(:)) max(a(:))]);
        colorbar;
        set(gca, 'ydir' , 'reverse'); % so that up is up, down is down
        title(gca,['t = ',num2str((ii*timeStep-tOffset1)/1000),' s']);
        set(gca, 'xtick', [] );
        set(gca, 'ytick', [] );
        box off;
        % errobar of azimuth
        axes('unit','pixels','pos',[60+100*ii+500*(ii-1) 200 315 100]);
        y_azimuth_mean = mean(spk_data_bin_mean_rate_trans{k}(:,:,ii),1);
        y_azimuth_std =std(spk_data_bin_mean_rate_trans{k}(:,:,ii),1);
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
        y_elevation_mean=mean(spk_data_bin_mean_rate_trans{k}(:,:,ii),2);
        y_elevation_std =std(spk_data_bin_mean_rate_trans{k}(:,:,ii),0,2);
        y_elevation_ste =y_elevation_std/ sqrt(length(find( (temp_elevation==unique_elevation(j))&(temp_stimType==unique_stimType(k)) )) );
        errorbar(yEle,y_elevation_mean,y_elevation_ste,'k.-','markerSize',20,'linewidth',2);
        xlim([-1 1]);
        view(90,270);
        xlabel('Elevation');
        ylim([min(y_elevation_mean)-max(y_elevation_ste), max(y_elevation_mean)+max(y_elevation_ste)+2]);
        axis off;
        %         pause(0.1);
        
        % text on the figure
        axes('unit','pixels','pos',[60 730 1800 80]);
        xlim([0,100]);
        ylim([0,10]);
        text(15,0,'Color contour maps','fontsize',28);
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
        text(11,0,num2str(maxSpkRealMean(k)),'fontsize',15);
        text(0,3,['Spon max FR(Hz): ',num2str(maxSpkSponMean)],'fontsize',15);
        text(0,6,['Spon mean FR(Hz): ',num2str(meanSpkSponMean)],'fontsize',15);
        axis off;
        
        drawnow;
        
        
        F=getframe(gcf);
        I=frame2im(F);
        [I,map]=rgb2ind(I,256);
        
        if Protocol == DIRECTION_TUNING_3D
            ss = [str, '_T'];
            file = ['Z:\LBY\Recording data\Qiaoqiao\3D_Tuning\Translation\' ss,'.gif'];
        elseif Protocol == ROTATION_TUNING_3D
            ss = [str, '_R'];
            file = ['Z:\LBY\Recording data\Qiaoqiao\3D_Tuning\Translation\' ss,'.gif'];
        end
        
        if pic_num == 1
            imwrite(I,map, file,'gif', 'Loopcount',inf,'DelayTime',0.1);
        else
            imwrite(I,map, file,'gif','WriteMode','append','DelayTime',0.1);
        end
        pic_num = pic_num + 1;
    end
    
    
    
end