%% First select proper directories and cells for exclusion
% datadir = 'Z:\Data\Tempo\Baskin\Analysis\CDD_split_behavior\LongDurCells\'
% datadir = 'Z:\Data\Tempo\Baskin\Analysis\CDD_split_behavior\ShortDurCells\'
datadir = 'Z:\Data\Tempo\Baskin\Analysis\CDD_split_behavior\ShortDurDelayCells\'
d = ls(datadir);
d = d(3:end,:);
% exclude_cells = [1 4 9 15 17 31 40 43 47 51 52 53 57 60]; %for long dur cells
% exclude_cells = [0]; %for shortdur cells
% exclude_cells = [2 15 18]; %for Baskin shortdur+delay cells
bask_longdur_exclude = {'m4c243r6','m4c250r5','m4c257r6','m4c262r5',...
    'm4c266r5','m4c267r7','m4c276r5','m4c284r5','m4c285r5',...
    'm4c186r5','m4c193r5','m4c204r6','m4c215r7','m4c221r6',... %
    'm4c239r5','m4c252r5','m4c253r5','m4c259r5','m4c267r7','m4c276r5','m4c282r5',...
    'm4c283r5','m4c290r6','m4c293r5'};
bask_shortdur_exclude = {};
bask_shortdur_delay_exclude = {'m4c313r5','m4c356r5','m4c361r5','m4c390r5','m4c401r5','m4c412r5','m4c413r5','m4c414r5','m4c423r5','m4c425r5','m4c426r5','m4c428r5'};
data_exclude = bask_shortdur_delay_exclude; %set this to the monkey and task version whose data you're combining

%% Now load data
coher = 2*[0 2 4 8 16];
cum_lort = []; cum_hirt = []; cum_loisi = []; cum_hiisi = [];
for i = 1:size(d,1)
    if ~ismember(d(i,4:11), data_exclude)
        c = importdata(sprintf('%s%s',datadir,d(i,:)));
        cum_lort(:,:,end+1) = c(2:4,:); %indices of mean rt at each coherence at each of 4+1(avg) conditions
        cum_hirt(:,:,end+1) = c(6:8,:);
        cum_loisi(:,:,end+1) = c(10:12,:); 
        cum_hiisi(:,:,end+1) = c(14:16,:); 
    end
end
cum_lort = cum_lort(:,:,2:end); %not sure why it's making a matrix of zeros
cum_hirt = cum_hirt(:,:,2:end); %not sure why it's making a matrix of zeros
cum_loisi = cum_loisi(:,:,2:end); %not sure why it's making a matrix of zeros
cum_hiisi = cum_hiisi(:,:,2:end); %not sure why it's making a matrix of zeros
n_cells = size(cum_lort,3);
LineTypes = {'b-','r-','g-';'b--','r--','g--'};
LineSymbols = {'bo','r*','g>'};
LineColors = {'b','r','g'};

handl(1) = figure; hold on;%% plot them all on one graph
set(handl(1),'PaperPosition', [.2 .2 8 10.7], 'Position', [250 50 500 573], 'Name', 'Group behavior: Split by RT or ISI');

%% first plot RT-split behavior
subplot(211); hold on;
mean_lort = mean(cum_lort,3);
se_lort = std(cum_lort,0,3)./sqrt(n_cells);
mean_hirt = mean(cum_hirt,3);
se_hirt = std(cum_hirt,0,3)./sqrt(n_cells);
for k = 1:3
    %fit lort and hirt
    fit_data(:,1) = coher';
    fit_data(:,2) = mean_lort(k,:);
    fit_data(:,3) = size(cum_lort,3);
    fit_x = [0:0.1:max(coher)];
    [group_alpha_lort(k) group_beta_lort(k) group_gamma_lort(k)] = weibull_bs_fit(fit_data);
%     [group_alpha_lort(k) group_beta_lort(k)] = logistic_fit(fit_data);
    fit_lort(k,:) = weibull_bs_curve(fit_x,[group_alpha_lort(k) group_beta_lort(k) group_gamma_lort(k)]);
%     fit_lort(k,:) = logistic_curve(fit_lort,[group_alpha_lort(k) group_beta_lort(k)]);
    fit_data(:,1) = coher';
    fit_data(:,2) = mean_hirt(k,:);
    fit_data(:,3) = size(cum_hirt,3);
    fit_x = [0:0.1:max(coher)];
    [group_alpha_hirt(k) group_beta_hirt(k) group_gamma_hirt(k)] = weibull_bs_fit(fit_data);
%     [group_alpha_hirt(k) group_beta_hirt(k)] = logistic_fit(fit_data);
    fit_hirt(k,:) = weibull_bs_curve(fit_x,[group_alpha_hirt(k) group_beta_hirt(k) group_gamma_hirt(k)]);
%     fit_hirt(k,:) = logistic_curve(fit_hirt,[group_alpha_hirt(k) group_beta_hirt(k)]);
    %now plot
    legh(k) = errorbar(coher, mean_lort(k,:), se_lort(k,:), LineSymbols{k});
    plot(fit_x, fit_lort(k,:), LineColors{1,k});
    errorbar(coher, mean_hirt(k,:), se_hirt(k,:), LineSymbols{k});
    plot(fit_x, fit_hirt(k,:), LineTypes{2,k});
end
legh(4)=plot([0 0],[-1 -1],'k-'); legh(5)=plot([0 0],[-1 -1],'k--');
xlim([-2 34]); ylim([0 1]); 
legend(legh,{'Invalid','Neutral','Valid','Short RT','Long RT'},'Location','SouthEast','Box','Off');
title('Behav split by RT')

%% now plot isi-split behavior
subplot(212); hold on;
mean_loisi = mean(cum_loisi,3);
se_loisi = std(cum_loisi,0,3)./sqrt(n_cells);
mean_hiisi = mean(cum_hiisi,3);
se_hiisi = std(cum_hiisi,0,3)./sqrt(n_cells);
for k = 1:3
    %fit loisi and hiisi
    fit_data(:,1) = coher';
    fit_data(:,2) = mean_loisi(k,:);
    fit_data(:,3) = size(cum_loisi,3);
    fit_x = [0:0.1:max(coher)];
    [group_alpha_loisi(k) group_beta_loisi(k) group_gamma_loisi(k)] = weibull_bs_fit(fit_data);
%     [group_alpha_loisi(k) group_beta_loisi(k)] = logistic_fit(fit_data);
    fit_loisi(k,:) = weibull_bs_curve(fit_x,[group_alpha_loisi(k) group_beta_loisi(k) group_gamma_loisi(k)]);
%     fit_loisi(k,:) = logistic_curve(fit_loisi,[group_alpha_loisi(k) group_beta_loisi(k)]);
    fit_data(:,1) = coher';
    fit_data(:,2) = mean_hiisi(k,:);
    fit_data(:,3) = size(cum_hiisi,3);
    fit_x = [0:0.1:max(coher)];
    [group_alpha_hiisi(k) group_beta_hiisi(k) group_gamma_hiisi(k)] = weibull_bs_fit(fit_data);
%     [group_alpha_hiisi(k) group_beta_hiisi(k)] = logistic_fit(fit_data);
    fit_hiisi(k,:) = weibull_bs_curve(fit_x,[group_alpha_hiisi(k) group_beta_hiisi(k) group_gamma_hiisi(k)]);
%     fit_hiisi(k,:) = logistic_curve(fit_hiisi,[group_alpha_hiisi(k) group_beta_hiisi(k)]);
    %now plot
    legh(k) = errorbar(coher, mean_loisi(k,:), se_loisi(k,:), LineSymbols{k});
    plot(fit_x, fit_loisi(k,:), LineColors{1,k});
    errorbar(coher, mean_hiisi(k,:), se_hiisi(k,:), LineSymbols{k});
    plot(fit_x, fit_hiisi(k,:), LineTypes{2,k});

end
legh(4)=plot([0 0],[-1 -1],'k-'); legh(5)=plot([0 0],[-1 -1],'k--');
xlim([-2 34]); ylim([0 1]); 
legend(legh,{'Invalid','Neutral','Valid','Short ISI','Long ISI'},'Location','SouthEast','Box','Off');
title('Behav split by ISI');

%% save data to text file
temp = [coher' mean_lort' mean_hirt'; coher' se_lort' se_hirt'; fit_x' fit_lort' fit_hirt'];
save('Z:\LabTools\Matlab\TEMPO_Analysis\ProtocolSpecific\CuedDirectionDiscrim\group_psych_splitbyRT_shortdur_delay.txt', '-ascii', 'temp')
temp = [coher' mean_loisi' mean_hiisi'; coher' se_loisi' se_hiisi'; fit_x' fit_loisi' fit_hiisi'];
save('Z:\LabTools\Matlab\TEMPO_Analysis\ProtocolSpecific\CuedDirectionDiscrim\group_psych_splitbyISI_shortdur_delay.txt', '-ascii', 'temp')
