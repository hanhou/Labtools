spu_exclude = [13 21 22 24]; %runs that don't use coher [0 2 4 8 16]
bask_exclude = [34 37 41 45 46 47 51 54 55];
bask_longdur_exclude = {'m4c243r6','m4c250r5','m4c257r6','m4c262r5',...
    'm4c266r5','m4c267r7','m4c276r5','m4c284r5','m4c285r5',...
    'm4c186r5','m4c193r5','m4c204r6','m4c215r7','m4c221r6',... %
    'm4c239r5','m4c252r5','m4c253r5','m4c259r5','m4c267r7','m4c276r5','m4c282r5',...
    'm4c283r5','m4c290r6','m4c293r5'};
bask_shortdur_exclude = {};
bask_shortdur_delay_exclude = {'m4c313r5','m4c356r5','m4c361r5','m4c390r5','m4c401r5','m4c412r5','m4c413r5','m4c414r5','m4c423r5','m4c425r5','m4c426r5','m4c428r5'};
data_exclude = bask_shortdur_delay_exclude; %set this to the monkey and task version whose data you're combining
% datadir = 'Z:\Data\Tempo\Spumoni\Analysis\1D-Psychophysics\';
% datadir = 'Z:\Data\Tempo\Baskin\Analysis\1D-Psychophysics\LongDurCells\';
% datadir = 'Z:\Data\Tempo\Baskin\Analysis\1D-Psychophysics\ShortDurCells\';
datadir = 'Z:\Data\Tempo\Baskin\Analysis\1D-Psychophysics\ShortDurDelayCells\';
temp = ls(datadir);
temp = temp(3:end,:);
d = [];
for i = 1:size(temp,1)
%     include = 1;
%     for j = 1:length(data_exclude)
%         if regexp(temp(i,:), sprintf('r%d',data_exclude(j))) %if that run should be excluded
%             include = 0;
%         end
%     end
%     if include, d = [d; temp(i,:)]; end    
%     if ~ismember(i, data_exclude), d = [d; temp(i,:)]; end
    if ~ismember(temp(i,7:14), data_exclude)
        d = [d; temp(i,:)]; 
    end
end
coher = 2*[0 2 4 8 16];
cum_psychdata = [];
for i = 1:size(d,1)
    c = importdata(sprintf('%s%s',datadir,d(i,:)));
    cum_psychdata(:,:,end+1) = c(2:4,1:5); %indices of pct correct at each validity
end
plot_colors = {'bs','ro','g>'};
plot_lines = {'b','r','g'};
plot_linecolor = {'bs-','ro-','g>-'};
mean_psychdata = mean(cum_psychdata,3);
std_psychdata = std(cum_psychdata,0,3);
figure; hold on;
clear fit_y;
for i = 1:3
    fit_data(:,1) = coher';
    fit_data(:,2) = mean_psychdata(i,:);
    fit_data(:,3) = size(cum_psychdata,3);
    fit_x = [0:0.1:max(coher)];
    [group_alpha(i) group_beta(i) group_gamma(i)] = weibull_bs_fit(fit_data);
%     [group_alpha(i) group_beta(i)] = logistic_fit(fit_data);
    fit_y(i,:) = weibull_bs_curve(fit_x,[group_alpha(i) group_beta(i) group_gamma(i)]);
%     fit_y(i,:) = logistic_curve(fit_x,[group_alpha(i) group_beta(i)]);
    errorbar(coher,mean_psychdata(i,:),std_psychdata(i,:)./sqrt(size(d,1)),plot_colors{i});
    plot(fit_x,fit_y(i,:),plot_lines{i});
end
axis('tight'); ylim([0.2 1]); axis('manual');
plot(xlim,[0.5 0.5],'k:')
for i = 1:3
    h(i) = plot(-1,-1,plot_linecolor{i});
end
l = legend(h,'Invalid','Neutral','Valid','Location','SouthEast');
set(l,'box','off')
xlabel('Coherence');
ylabel('Percent Correct');
title(sprintf('(n=%d)',size(d,1)));

temp = [coher' mean_psychdata'; coher' (std_psychdata./sqrt(size(d,1)))'; fit_x' fit_y'];
save('Z:\LabTools\Matlab\TEMPO_Analysis\ProtocolSpecific\CuedDirectionDiscrim\group_psych_shortdur_delay.txt', '-ascii', 'temp')
