figure;
xlabels={'RT1','RT2','Holding Time','HoldingTime','Reward'};
maxJ = 5;
maxI = 7;
for i=1:maxI
    for j=1:maxJ
        subplot(maxI,maxJ,maxJ*(i-1)+j);
        plot(1,1);
        if i==maxI
            xlabel(xlabels{j});
        end
    end
end