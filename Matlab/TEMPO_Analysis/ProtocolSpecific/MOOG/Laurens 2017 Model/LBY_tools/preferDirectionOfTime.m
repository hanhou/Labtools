% plot figures of preferred direction as a function of time
% LBY 20170325

% ------ fig.60 plot preferred direction as time ------%

for k = 1:length(unique_stimType)
    figure(60+k);
    set(gcf,'pos',[0 0 1900 1000]);
    clf;
    subplot(1,2,1);
%     plot(1:nBins(1,1),preDirAcrossTBin{k}(:,1)
    plot(preDirAcrossTBin{k}(:,1),preDirAcrossTBin{k}(:,2),'r-o');
    
    subplot(1,2,2);
    
    
end