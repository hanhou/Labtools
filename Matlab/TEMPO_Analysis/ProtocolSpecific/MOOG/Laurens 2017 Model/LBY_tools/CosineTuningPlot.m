% plot figures of Firing Rates (FR) as a function of angles between each
% directions & preferred direction across time bins
% LBY 20170324

global PSTH;

inx_row = floor(sqrt(nBins(1,1)/2));
inx_col = round(nBins(1,1)/inx_row);
% ------ fig.60 plot cosine tuning of directions across bins ------%

for k = 1:length(unique_stimType)
    figure(60+k);
    set(gcf,'pos',[0 0 1900 1000]);
    clf;
    index = 1;
    for nn = 1:nBins(1,1)
        subplot(inx_row,inx_col,index);
        for j = 1:length(unique_elevation)
            for i = 1:length(unique_azimuth)
%                 plot(sin(angleDiffs{k}(j,i,nn)),squeeze(PSTH.spk_data_bin_mean_rate{k}(j,i,nn)),'r.');
                 plot(angleDiffs{k}(j,i,nn),squeeze(PSTH.spk_data_bin_mean_rate{k}(j,i,nn)),'r.');
                hold on;
            end
        end
        set(gca,'xlim',[0 180]);
        axis off;
        index = index + 1;
    end
    
    %text on the figure
    axes('unit','pixels','pos',[60 850 1800 80]);
    xlim([0,100]);
    ylim([0,10]);
    if Protocol == DIRECTION_TUNING_3D
        text(30,10,'Cosine tuning (T)','fontsize',20);
    elseif Protocol == ROTATION_TUNING_3D
        text(30,10,'Cosine tuning(R)','fontsize',20);
    end
 
    FileNameTemp = num2str(FILE);
    FileNameTemp =  FileNameTemp(1:end);
    str = [FileNameTemp,'_Ch' num2str(SpikeChan)];
    str1 = [FileNameTemp,'\_Ch' num2str(SpikeChan),'    ',stimType{k}];
    text(70,0,str1,'fontsize',18);
    axis off;
    
% save the figure
    str6 = [str,'_cosine_tuning_',stimType{k}];
    set(gcf,'paperpositionmode','auto');
    switch Protocol
        case DIRECTION_TUNING_3D
            ss = [str6, '_T'];
            saveas(60+k,['Z:\LBY\Recording data\Qiaoqiao\3D_Tuning\Translation\' ss], 'emf');
        case ROTATION_TUNING_3D
            ss = [str6, '_R'];
            saveas(60+k,['Z:\LBY\Recording data\Qiaoqiao\3D_Tuning\Rotation\' ss], 'emf');
    end
end