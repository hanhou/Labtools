LONGDUR = 1;
SHORTDUR = 2;
SHORTDURDELAY = 3;

% % for Baskin's shortdur data use the following
% analysis_dir = 'Z:\Data\Tempo\Baskin\Analysis\LIP_PSTH\ShortDur';
% d = ls(analysis_dir);
% d = d(3:end,:); %removes '.' and '..' from top of list
% task = SHORTDUR;

%for Baskin's shortdur+delay data use the following
analysis_dir = 'Z:\Data\Tempo\Baskin\Analysis\LIP_PSTH\ShortDurDelay';
d = ls(analysis_dir);
d = d(3:end,:); %removes '.' and '..' from top of list
% bask_shortdur_delay_exclude = {'m4c313r5','m4c356r5','m4c361r5','m4c390r5','m4c401r5','m4c412r5','m4c413r5','m4c414r5','m4c423r5','m4c425r5','m4c426r5','m4c428r5', ...
%     'm4c439r5','m4c441r5','m4c444r5','m4c446r5','m4c447r5','m4c448r5','m4c451r5','m4c457r5'};
bask_shortdur_delay_exclude = {'m4c312r5', 'm4c317r5', 'm4c356r5', 'm4c383r5', 'm4c388r8', 'm4c412r5', ...
    'm4c413r5', 'm4c414r5', 'm4c428r5', 'm4c434r5', 'm4c439r5', 'm4c441r5', 'm4c444r5', 'm4c446r5', 'm4c447r5',...
    'm4c448r5', 'm4c451r5', 'm4c452r5',  'm4c457r5'};  

data_exclude = bask_shortdur_delay_exclude;
yellowbar_start = 0;
yellowbar_end = 250;
greybar_start = 250;
greybar_end = 500;
task = SHORTDURDELAY;

% % %for Baskin's full dur(1000ms) data use the following:
% analysis_dir = 'Z:\Data\Tempo\Baskin\Analysis\LIP_PSTH\LongDur';
% d = ls(analysis_dir);
% d = d(3:end,:); %removes '.' and '..' 
% % d = d([1:60 62:length(d)],:) %KLUGE to remove cell 61 which has problems 
% data_exclude = {'m4c243r6','m4c250r5','m4c257r6','m4c262r5',...
%     'm4c266r5','m4c267r7','m4c276r5','m4c284r5','m4c285r5',...
%     'm4c186r5','m4c193r5','m4c204r6','m4c215r7','m4c221r6',... 
%     'm4c239r5','m4c252r5','m4c253r5','m4c259r5','m4c267r7','m4c276r5','m4c282r5',...
%     'm4c283r5','m4c290r6','m4c293r5'};
% yellowbar_start = 0;
% yellowbar_end = 1000;
% greybar_start = 1000;
% greybar_end = 1000;
% task = LONGDUR;

%% set up matrices and load data
reload = 1; 

if reload
cum_file = []; cum_coher = []; cum_postcue_combined_psth = repmat({[]},[1,3]); cum_normval = []; cum_num_trials = [];
cum_postcue_psth = repmat({[]},[2,3]); cum_prestim_psth = repmat({[]},[5,2,3]); cum_sacc_psth = repmat({[]},[5,2,3]); 
cum_co_postcue_psth = repmat({[]},[2,3]); cum_co_prestim_psth = repmat({[]},[2,3]); cum_co_sacc_psth = repmat({[]},[2,3]); 
cum_ungrouped_prestim_psth = repmat({[]},[5,2,3,2]); cum_ungrouped_sacc_psth = repmat({[]},[5,2,3,2]); 
cum_lort_prestim_psth = repmat({[]},[5,2,3]); cum_hirt_prestim_psth = repmat({[]},[5,2,3]); 
cum_loisi_prestim_psth = repmat({[]},[5,2,3]); cum_hiisi_prestim_psth = repmat({[]},[5,2,3]); 
cum_zero_pct_xing = []; cum_co_postcue_combined_psth = repmat({[]},[1,3]); cum_motgrouped_prestim_psth = repmat({[]},[5,3,2]); 
% for m = 1:length(d)
%     load(sprintf('%s\\%s',analysis_dir,d(m,:)));
%     cum_file{length(cum_file)+1} = file;
%     cum_coher = [cum_coher; coher];
%     cum_normval = [cum_normval; normval];
%     for k = 1:3 %length(unique_cue_dir_type)
%         cum_postcue_combined_psth{k}(end+1,:) = sm_postcue_combined_psth{k};
%         for j = 1:2 %prefchoice = 1, nullchoice = 2
%             cum_postcue_psth{j,k}(end+1,:) = sm_postcue_psth{j,k};
%             cum_co_postcue_psth{j,k}(end+1,:) = sm_co_postcue_psth{j,k};
%             cum_co_prestim_psth{j,k}(end+1,:) = sm_co_prestim_psth{j,k};
%             cum_co_sacc_psth{j,k}(end+1,:) = sm_co_sacc_psth{j,k};
%             for i = 1:length(coher)
%                 cum_prestim_psth{i,j,k}(end+1,:) = sm_prestim_psth{i,j,k};
%                 cum_sacc_psth{i,j,k}(end+1,:) = sm_sacc_psth{i,j,k};
%             end
%         end
%     end
% end
%    
% SAVEFILE = 'Z:\LabTools\Matlab\TEMPO_Analysis\ProtocolSpecific\CuedDirectionDiscrim\cum_lip_psth_longdur_summary.mat';
% save(SAVEFILE, 'cum_file','cum_coher','cum_normval',...
%     'cum_postcue_combined_psth','cum_postcue_psth','cum_prestim_psth','cum_sacc_psth',...
%     'cum_co_postcue_psth','cum_co_prestim_psth','cum_co_sacc_psth');

unique_direction = [1 0]; %prefdir nulldir
stddev = 30;
buff = 3*stddev;
gaussfilt = normpdf([1:2*buff+1],buff+1,stddev); %gaussian filter 3 std.dev's wide
eff_width = sqrt(15^2+stddev^2)
lb = 100; %short for long_buffer... 100ms on either side
excluded_cells_ungrouped_data = cell(5,2,3,2);
for m = 1:size(d,1)
    load(sprintf('%s\\%s',analysis_dir,d(m,:)));
    %HAACK
    long_lort_sm_prestim_psth = long_sm_prestim_psth;
    long_hirt_sm_prestim_psth = long_sm_prestim_psth;
    long_loisi_sm_prestim_psth = long_sm_prestim_psth;
    long_hiisi_sm_prestim_psth = long_sm_prestim_psth;
    %HAACK
    cum_file{length(cum_file)+1} = file;
    cum_coher = [cum_coher; coher];
    cum_normval = [cum_normval; normval];
    cum_zero_pct_xing = [cum_zero_pct_xing; offset];
    cum_num_trials(m,:,:,:,:) = num_trials;
    for k = 1:3 %length(unique_cue_dir_type)
        temp = conv(gaussfilt, long_sm_postcue_combined_psth{k});
        cum_postcue_combined_psth{k}(end+1,:) = temp(lb+buff+1:end-lb-buff);
        if isequal(size(long_sm_co_postcue_psth{1,k}),size(long_sm_co_postcue_psth{2,k}))
            temp = mean([long_sm_co_postcue_psth{1,k};long_sm_co_postcue_psth{2,k}],1); %make a temporary quick avg between the psths in both directions
        else
            if size(long_sm_co_postcue_psth{1,k},1)==1 
                temp = mean([long_sm_co_postcue_psth{1,k};long_sm_co_postcue_psth{2,k}'],1); %make a temporary quick avg between the psths in both directions
            else
                temp = mean([long_sm_co_postcue_psth{1,k}';long_sm_co_postcue_psth{2,k}],1); %make a temporary quick avg between the psths in both directions
            end
        end
        temp = conv(gaussfilt, temp);
        cum_co_postcue_combined_psth{k}(end+1,:) = temp(lb+buff+1:end-lb-buff);
        for j = 1:2 %prefchoice = 1, nullchoice = 2
            temp = conv(gaussfilt, long_sm_postcue_combined_psth{k});
            cum_postcue_psth{j,k}(end+1,:) = temp(lb+buff+1:end-lb-buff);
            if k ~= 2
                temp = conv(gaussfilt, long_sm_co_postcue_psth{j,k});
                cum_co_postcue_psth{j,k}(end+1,:) = temp(lb+buff+1:end-lb-buff);
                temp = conv(gaussfilt, long_sm_co_prestim_psth{j,k});
                cum_co_prestim_psth{j,k}(end+1,:) = temp(lb+buff+1:end-lb-buff);
                temp = conv(gaussfilt, long_sm_co_sacc_psth{j,k});
                cum_co_sacc_psth{j,k}(end+1,:) = temp(lb+buff+1:end-lb-buff);
            end
            for i = 1:length(coher)
                temp = conv(gaussfilt, long_sm_prestim_psth{i,j,k});
                cum_prestim_psth{i,j,k}(end+1,:) = temp(lb+buff+1:end-lb-buff);
                temp = conv(gaussfilt, long_sm_sacc_psth{i,j,k});
                cum_sacc_psth{i,j,k}(end+1,:) = temp(lb+buff+1:end-lb-buff);
                temp = conv(gaussfilt, long_lort_sm_prestim_psth{i,j,k});
                cum_lort_prestim_psth{i,j,k}(end+1,:) = temp(lb+buff+1:end-lb-buff);
                temp = conv(gaussfilt, long_hirt_sm_prestim_psth{i,j,k});
                cum_hirt_prestim_psth{i,j,k}(end+1,:) = temp(lb+buff+1:end-lb-buff);
                temp = conv(gaussfilt, long_loisi_sm_prestim_psth{i,j,k});
                cum_loisi_prestim_psth{i,j,k}(end+1,:) = temp(lb+buff+1:end-lb-buff);
                temp = conv(gaussfilt, long_hiisi_sm_prestim_psth{i,j,k});
                cum_hiisi_prestim_psth{i,j,k}(end+1,:) = temp(lb+buff+1:end-lb-buff);
                for g = 1:2 %length(unique_direction)
                    temp = conv(gaussfilt, long_sm_ungrouped_prestim_psth{i,j,k,g});
                    if size(temp,1)>1, temp = temp'; end
                    cum_ungrouped_prestim_psth{i,j,k,g}(m,:) = temp(lb+buff+1:end-lb-buff);
                    temp = conv(gaussfilt, long_sm_ungrouped_sacc_psth{i,j,k,g});
                    if size(temp,1)>1, temp = temp'; end
                    cum_ungrouped_sacc_psth{i,j,k,g}(m,:) = temp(lb+buff+1:end-lb-buff);
                    if ~num_trials(i,j,k,g)
                        excluded_cells_ungrouped_data{i,j,k,g}(end+1) = m;
                    end
                    if j == 1 %do weighted avg of choice trials to get psths sorted by motion (and not choice)
                        temp = ( num_trials(i,1,k,g).*long_sm_ungrouped_prestim_psth{i,1,k,g} + ...
                                 num_trials(i,2,k,g).*long_sm_ungrouped_prestim_psth{i,2,k,g} ) ...
                               ./ (num_trials(i,1,k,g)+num_trials(i,2,k,g)); %avg of psth for each choice, weighted by num of choices in each dir
                        temp = conv(gaussfilt, temp);
                        if size(temp,1)>1, temp = temp'; end
                        cum_motgrouped_prestim_psth{i,k,g}(m,:) = temp(lb+buff+1:end-lb-buff);
                    end
                end
            end
        end
    end
end


end


% DATAFILE = 'Z:\LabTools\Matlab\TEMPO_Analysis\ProtocolSpecific\CuedDirectionDiscrim\cum_lip_psth_TEST.mat'

%% define constants
% load(DATAFILE);
n_total = length(cum_file)
gooddata = true(n_total,1); %these indicies indicate the subset of datasets that are well-isolated and show decent behavior
for i = 1:n_total
    if ismember(cum_file{i}(1:end-4), data_exclude)
        gooddata(i) = 0;
    end
end
gooddata_co = gooddata;
for i = 1:length(gooddata_co) %exclude cells where the monkey didn't err at least once on both directions on the cueonly trials
    if ( isnan(cum_co_postcue_psth{1,1}(i,1)) || isnan(cum_co_postcue_psth{2,1}(i,1)) || ...
            isnan(cum_co_postcue_psth{2,3}(i,1)) || isnan(cum_co_postcue_psth{2,3}(i,1)) ) 
        gooddata_co(i) = 0;
    end
end
gooddata_co = gooddata; %maybe the above criterion is TOO restrictive... won't let you do paired t-tests though
n = sum(gooddata);
n_co = sum(gooddata_co);

coh_lines = {'r','m','b','c','g'};
coh_dashedlines = {'r--','m--','b--','c--','g--'};
coh_dottedlines = {'r:','m:','b:','c:','g:'};

unique_cue_dir_type = [-1 0 1]; %-1 = NullCue, 0 = Neutral Cue, 1 = PrefCue
unique_cue_dir_type_names = {'NullCue','Neutral','PrefCue'};
unique_coherence = cum_coher(1,:);
unique_direction = [1 2]; %1=PrefDir motion, 2=NullDir motion
%these timings are determined in PSTH_CuedDirec
precue = 200; %time to show before the cue starts
postcue = 400; %time to show after cue starts
prestim = 300; %time to display before the visual stimulus
poststim = 1100; %time to display after the visual stimulus
presacc = 700; %time to show before the saccade starts
postsacc = 200; %time to show after saccade starts
precue_co = 200;
postcue_co = 200+150+1000;
prestim_co = -800; %200 before FP off
poststim_co = 1400;
postcue_x = [-precue:(size(cum_postcue_psth{1,1,1},2)-precue-1)];
prestim_x = [-prestim:(size(cum_prestim_psth{1,1,1},2)-prestim-1)];
postcue_co_x = [-precue_co:(length(cum_co_postcue_psth{1,1})-precue_co-1)];
prestim_co_x = [-200:(length(cum_co_prestim_psth{1,1})-200-1)]; %relative to FP offset

sacc_x = [-presacc:(size(cum_sacc_psth{1,1},2)-presacc-1)];
linetypes = {'b-','r-','g-'};
linetypesT = {'b-','r-','g-';'b:','r:','g:'};
linesymbols = {'bo-','r*-','g<-'};
linesymbolsT = {'bo-','r*-','g<-';'bo:','r*:','g<:'};
clear temp cued_sacc_onset co_sacc_onset combined_sacc_onset;
cued_sacc_onset = cell(length(unique_coherence), length(unique_cue_dir_type));
co_sacc_onset = cell(length(unique_cue_dir_type));
combined_sacc_onset = cell(length(unique_cue_dir_type));
for k = 1:length(unique_cue_dir_type)
    for i = 1:length(unique_coherence)
        temp(:,i) = cum_sacc_psth{i,1,k}(gooddata,presacc+1); %use this to normalize non-cueonly sorted by coherence trials
        cued_sacc_onset{i,k} = temp(:,i);
    end
    co_sacc_onset{k} = cum_sacc_psth{1,k}(gooddata,presacc+1); %use this to normalize cueonly trials
    combined_sacc_onset{k} = mean(temp,2); %this takes a mean of the saccade onsets across all coherences
    %use to normalize cue onset psths
end
normval = cum_normval(:,2); %use peak around target onset to normalize

%% for each cell, plot the saccadic response, averaged across coherence and cue type
figure;
cellsaccpsth = zeros(n_total,2,length(sacc_x));
for cellcounter = 1:size(cum_sacc_psth{1},1)
    cellsaccpsth(cellcounter,1,:) = nanmean([cum_sacc_psth{1,1,1}(cellcounter,:); cum_sacc_psth{2,1,1}(cellcounter,:); cum_sacc_psth{3,1,1}(cellcounter,:); cum_sacc_psth{4,1,1}(cellcounter,:); cum_sacc_psth{5,1,1}(cellcounter,:); ...
        cum_sacc_psth{1,1,2}(cellcounter,:); cum_sacc_psth{2,1,2}(cellcounter,:); cum_sacc_psth{3,1,2}(cellcounter,:); cum_sacc_psth{4,1,2}(cellcounter,:); cum_sacc_psth{5,1,2}(cellcounter,:); ...
        cum_sacc_psth{1,1,3}(cellcounter,:); cum_sacc_psth{2,1,3}(cellcounter,:); cum_sacc_psth{3,1,3}(cellcounter,:); cum_sacc_psth{4,1,3}(cellcounter,:); cum_sacc_psth{5,1,3}(cellcounter,:)],1);
    cellsaccpsth(cellcounter,2,:) = nanmean([cum_sacc_psth{1,2,1}(cellcounter,:); cum_sacc_psth{2,2,1}(cellcounter,:); cum_sacc_psth{3,2,1}(cellcounter,:); cum_sacc_psth{4,2,1}(cellcounter,:); cum_sacc_psth{5,2,1}(cellcounter,:); ...
        cum_sacc_psth{1,2,2}(cellcounter,:); cum_sacc_psth{2,2,2}(cellcounter,:); cum_sacc_psth{3,2,2}(cellcounter,:); cum_sacc_psth{4,2,2}(cellcounter,:); cum_sacc_psth{5,2,2}(cellcounter,:); ...
        cum_sacc_psth{1,2,3}(cellcounter,:); cum_sacc_psth{2,2,3}(cellcounter,:); cum_sacc_psth{3,2,3}(cellcounter,:); cum_sacc_psth{4,2,3}(cellcounter,:); cum_sacc_psth{5,2,3}(cellcounter,:)],1);
    subplot(5,ceil(n_total/5),cellcounter); hold on;
    plot(sacc_x, squeeze(cellsaccpsth(cellcounter,1,:)),'k-');
    plot(sacc_x, squeeze(cellsaccpsth(cellcounter,2,:)),'k--');
    axis tight;
    line([0 0], ylim, 'Color','r');
    title(cum_file{cellcounter});
end


%% first compute means and normalize means for all psth.  Normalization is
%to the target onset of the corresponding trial type
%also for those data sets where the monkey didn't make both possible
%choices for the cue only trials.  

%first preallocate matrices
mean_postcue_combined_psth = cell(length(unique_cue_dir_type),1);
mean_co_postcue_combined_psth = cell(length(unique_cue_dir_type),1);
mean_postcue_psth = cell(2,length(unique_cue_dir_type));
mean_co_postcue_psth = cell(2,length(unique_cue_dir_type));
mean_co_prestim_psth = cell(2,length(unique_cue_dir_type));
mean_co_sacc_psth = cell(2,length(unique_cue_dir_type));
mean_sacc_psth = cell(length(unique_coherence),2,length(unique_cue_dir_type));
mean_prestim_psth = cell(length(unique_coherence),2,length(unique_cue_dir_type));
mean_loisi_prestim_psth = cell(length(unique_coherence),2,length(unique_cue_dir_type));
mean_hiisi_prestim_psth = cell(length(unique_coherence),2,length(unique_cue_dir_type));
mean_lort_prestim_psth = cell(length(unique_coherence),2,length(unique_cue_dir_type));
mean_hirt_prestim_psth = cell(length(unique_coherence),2,length(unique_cue_dir_type));
mean_ungrouped_prestim_psth = cell(length(unique_coherence),2,length(unique_cue_dir_type),2);
mean_ungrouped_sacc_psth = cell(length(unique_coherence),2,length(unique_cue_dir_type),2);
mean_motgrouped_prestim_psth = cell(length(unique_coherence),length(unique_cue_dir_type),2);

for k = 1:length(unique_cue_dir_type)
    
    select = ~isnan(cum_postcue_combined_psth{k}(:,1));
    mean_postcue_combined_psth{k} = mean(cum_postcue_combined_psth{k}(gooddata & select,:),1);
%     normval_postcue_combined_psth{k} = repmat(normval,1,size(cum_postcue_combined_psth{k},2));
%     norm_postcue_combined_psth{k} = cum_postcue_combined_psth{k}./normval_postcue_combined_psth{k};
%     normmean_postcue_combined_psth{k} = mean(norm_postcue_combined_psth{k}(gooddata & select,:),1);
    
    if k ~=2
        select = ~isnan(cum_co_postcue_combined_psth{k}(:,1));
        mean_co_postcue_combined_psth{k} = mean(cum_co_postcue_combined_psth{k}(gooddata_co & select,:),1);
%         normval_co_postcue_combined_psth{k} = repmat(normval,1,size(cum_co_postcue_combined_psth{k},2));
%         norm_co_postcue_combined_psth{k} = cum_co_postcue_combined_psth{k}./normval_co_postcue_combined_psth{k};
%         normmean_co_postcue_combined_psth{k} = mean(norm_co_postcue_combined_psth{k}(gooddata_co & select,:),1);
    end
    for j = 1:2
        select = ~isnan(cum_postcue_psth{j,k}(:,1));
        mean_postcue_psth{j,k} = mean(cum_postcue_psth{j,k}(gooddata & select,:),1);
%         normval_postcue_psth{j,k} = repmat(normval,1,size(cum_postcue_psth{j,k},2));
%         norm_postcue_psth{j,k} = cum_postcue_psth{j,k}./normval_postcue_psth{j,k};
%         normmean_postcue_psth{j,k} = mean(norm_postcue_psth{j,k}(gooddata & select,:),1);

        if k ~= 2
            select = ~isnan(cum_co_postcue_psth{j,k}(:,1));
            mean_co_postcue_psth{j,k} = mean(cum_co_postcue_psth{j,k}(gooddata_co & select,:),1);
%             normval_co_postcue_psth{j,k} = repmat(normval,1,size(cum_co_postcue_psth{j,k},2));
%             norm_co_postcue_psth{j,k} = cum_co_postcue_psth{j,k}./normval_co_postcue_psth{j,k};
%             normmean_co_postcue_psth{j,k} = mean(norm_co_postcue_psth{j,k}(gooddata_co & select,:),1);

            select = ~isnan(cum_co_prestim_psth{j,k}(:,1));
            mean_co_prestim_psth{j,k} = mean(cum_co_prestim_psth{j,k}(gooddata_co & select,:),1);
%             normval_co_prestim_psth{j,k} = repmat(normval,1,size(cum_co_prestim_psth{j,k},2));
%             norm_co_prestim_psth{j,k} = cum_co_prestim_psth{j,k}./normval_co_prestim_psth{j,k};
%             normmean_co_prestim_psth{j,k} = mean(norm_co_prestim_psth{j,k}(gooddata_co & select,:),1);

            select = ~isnan(cum_co_sacc_psth{j,k}(:,1));
            mean_co_sacc_psth{j,k} = mean(cum_co_sacc_psth{j,k}(gooddata_co & select,:),1);
%             normval_co_sacc_psth{j,k} = repmat(normval,1,size(cum_co_sacc_psth{j,k},2));
%             norm_co_sacc_psth{j,k} = cum_co_sacc_psth{j,k}./normval_co_sacc_psth{j,k};
%             normmean_co_sacc_psth{j,k} = mean(norm_co_sacc_psth{j,k}(gooddata_co & select,:),1);
        end

        for i = 1:length(unique_coherence)
            select = ~isnan(cum_prestim_psth{i,j,k}(:,1));
            mean_prestim_psth{i,j,k} = mean(cum_prestim_psth{i,j,k}(gooddata & select,:),1);
%             normval_prestim_psth{i,j,k} = repmat(normval,1,size(cum_prestim_psth{i,j,k},2));
%             norm_prestim_psth{i,j,k} = cum_prestim_psth{i,j,k}./normval_prestim_psth{i,j,k};
%             normmean_prestim_psth{i,j,k} = mean(norm_prestim_psth{i,j,k}(gooddata & select,:),1);

            select = ~isnan(cum_sacc_psth{i,j,k}(:,1));
            mean_sacc_psth{i,j,k} = mean(cum_sacc_psth{i,j,k}(gooddata & select,:),1);
%             normval_sacc_psth{i,j,k} = repmat(normval,1,size(cum_sacc_psth{i,j,k},2));
%             norm_sacc_psth{i,j,k} = cum_sacc_psth{i,j,k}./normval_sacc_psth{i,j,k};
%             normmean_sacc_psth{i,j,k} = mean(norm_sacc_psth{i,j,k}(gooddata & select,:),1);
            
            select = ~isnan(cum_lort_prestim_psth{i,j,k}(:,1));
            mean_lort_prestim_psth{i,j,k} = mean(cum_lort_prestim_psth{i,j,k}(gooddata & select,:),1);
%             normval_lort_prestim_psth{i,j,k} = repmat(normval,1,size(cum_lort_prestim_psth{i,j,k},2));
%             norm_lort_prestim_psth{i,j,k} = cum_lort_prestim_psth{i,j,k}./normval_lort_prestim_psth{i,j,k};
%             normmean_lort_prestim_psth{i,j,k} = mean(norm_lort_prestim_psth{i,j,k}(gooddata & select,:),1);
            
            select = ~isnan(cum_hirt_prestim_psth{i,j,k}(:,1));
            mean_hirt_prestim_psth{i,j,k} = mean(cum_hirt_prestim_psth{i,j,k}(gooddata & select,:),1);
%             normval_hirt_prestim_psth{i,j,k} = repmat(normval,1,size(cum_hirt_prestim_psth{i,j,k},2));
%             norm_hirt_prestim_psth{i,j,k} = cum_hirt_prestim_psth{i,j,k}./normval_hirt_prestim_psth{i,j,k};
%             normmean_hirt_prestim_psth{i,j,k} = mean(norm_hirt_prestim_psth{i,j,k}(gooddata & select,:),1);

            select = ~isnan(cum_loisi_prestim_psth{i,j,k}(:,1));
            mean_loisi_prestim_psth{i,j,k} = mean(cum_loisi_prestim_psth{i,j,k}(gooddata & select,:),1);
%             normval_loisi_prestim_psth{i,j,k} = repmat(normval,1,size(cum_loisi_prestim_psth{i,j,k},2));
%             norm_loisi_prestim_psth{i,j,k} = cum_loisi_prestim_psth{i,j,k}./normval_loisi_prestim_psth{i,j,k};
%             normmean_loisi_prestim_psth{i,j,k} = mean(norm_loisi_prestim_psth{i,j,k}(gooddata & select,:),1);
            
            select = ~isnan(cum_hiisi_prestim_psth{i,j,k}(:,1));
            mean_hiisi_prestim_psth{i,j,k} = mean(cum_hiisi_prestim_psth{i,j,k}(gooddata & select,:),1);
%             normval_hiisi_prestim_psth{i,j,k} = repmat(normval,1,size(cum_hiisi_prestim_psth{i,j,k},2));
%             norm_hiisi_prestim_psth{i,j,k} = cum_hiisi_prestim_psth{i,j,k}./normval_hiisi_prestim_psth{i,j,k};
%             normmean_hiisi_prestim_psth{i,j,k} = mean(norm_hiisi_prestim_psth{i,j,k}(gooddata & select,:),1);

            for g = 1:2
                temp = gooddata;
                temp(excluded_cells_ungrouped_data{i,j,k,g}) = 0; %exclude those cells that lacked the corresponding condition
                
                mean_ungrouped_prestim_psth{i,j,k,g} = mean(cum_ungrouped_prestim_psth{i,j,k,g}(temp,:),1);
%                 normval_ungrouped_prestim_psth{i,j,k,g} = repmat(normval,1,size(cum_ungrouped_prestim_psth{i,j,k,g},2));
%                 norm_ungrouped_prestim_psth{i,j,k,g} = cum_ungrouped_prestim_psth{i,j,k,g}./normval_ungrouped_prestim_psth{i,j,k,g};
%                 normmean_ungrouped_prestim_psth{i,j,k,g} = mean(norm_ungrouped_prestim_psth{i,j,k,g}(temp,:),1);

                mean_ungrouped_sacc_psth{i,j,k,g} = mean(cum_ungrouped_sacc_psth{i,j,k,g}(temp,:),1);
%                 normval_ungrouped_sacc_psth{i,j,k,g} = repmat(normval,1,size(cum_ungrouped_sacc_psth{i,j,k,g},2));
%                 norm_ungrouped_sacc_psth{i,j,k,g} = cum_ungrouped_sacc_psth{i,j,k,g}./normval_ungrouped_sacc_psth{i,j,k,g};
%                 normmean_ungrouped_sacc_psth{i,j,k,g} = mean(norm_ungrouped_sacc_psth{i,j,k,g}(temp,:),1);
                
                if j == 1
                    select = ~isnan(cum_motgrouped_prestim_psth{i,k,g}(:,1));
                    mean_motgrouped_prestim_psth{i,k,g} = mean(cum_motgrouped_prestim_psth{i,k,g}((select & gooddata),:),1);
                end                    
            end
        end
    end
end

%% try to assess when the motion integration starts: 
%for each cell, take the neutral cue and subtract off the T2 motion psth
%from the T1 motion psth, then average and compute the 95% CI (by bootstrap?) and see
%when that excludes 0.  Maybe do this separately by coherence and then
%average?  Only use correct trials so that the 
neutdiff_psth = {[],[],[],[],[]};
for m = 1:length(cum_file)
    if gooddata(m)
        for i = 1:5
            neutdiff_psth{i}(end+1,:) = cum_ungrouped_prestim_psth{i,1,2,1}(m,:)-cum_ungrouped_prestim_psth{i,2,2,2}(m,:); 
%             neutdiff_psth{i}(end+1,:) = cum_motgrouped_prestim_psth{i,2,1}(m,:)-cum_motgrouped_prestim_psth{i,2,2}(m,:); 
        end
    end
end
%now do bootstrap over the data to get a multiple-comparison corrected
%95%CI of the mean (ie, 99%CI given 5 coherences).
nboot = 2000;
bootmeandiff_psth = cell(length(coher),1);
meandiff_psth = zeros(length(coher),length(prestim_x));
sortbootmeandiff_psth = cell(length(coher),1);
CIhi_meandiff_psth = zeros(length(coher),length(prestim_x));
CIlo_meandiff_psth = zeros(length(coher),length(prestim_x));
diff_latency = zeros(length(coher),1);
for b = 1:nboot
    bootshuff = ceil(sum(gooddata).*rand(sum(gooddata),1));
    for i = 1:length(coher)
        bootmeandiff_psth{i}(b,:) = mean(neutdiff_psth{i}(bootshuff,:));
    end
end
figure; hold on
for i = 1:5
    meandiff_psth(i,:) = nanmean(neutdiff_psth{i}); %averages across cells
    sortbootmeandiff_psth{i} = sort(bootmeandiff_psth{i});
    CIhi_meandiff_psth(i,:) = sortbootmeandiff_psth{i}(ceil(0.995*nboot),:);
    CIlo_meandiff_psth(i,:) = sortbootmeandiff_psth{i}(floor(0.005*nboot),:);
    plot(prestim_x, meandiff_psth(i,:), coh_lines{i});
    plot(prestim_x, CIhi_meandiff_psth(i,:), coh_dottedlines{i});
    plot(prestim_x, CIlo_meandiff_psth(i,:), coh_dottedlines{i});
    diff_latency(i) = prestim_x(min([length(prestim_x) (find(CIlo_meandiff_psth(i,1:701)<=0,1,'last')+1)]));
end
plot(xlim,[0 0],'k-');
xlabel('Time post motion onset'); ylabel('Delta firing rate (Hz)');
diff_latency %this would print the latencies in order of coherence
startperiod = min(diff_latency);
endperiod = startperiod + yellowbar_end;

%% setup data(y) and covariate(x) matrices for the ramp and baserate portions
%of the cell-by-cell and population data.  This prepares 4 regressions.
begin_ramp = startperiod; end_ramp = endperiod;
x_popramp = []; y_popramp = []; y_poprampnorm = [];
x_cellramp = []; y_cellramp = []; cellrampavg = []; 
rampinds = find( (prestim_x >= begin_ramp) & (prestim_x <= end_ramp) );

begin_base = -300; end_base = 200;
x_popbase = []; y_popbase = []; y_popbasenorm = []; 
x_cellbase = []; y_cellbase = []; cellbaseavg = [];
baseinds = find( (prestim_x >= begin_base) & (prestim_x <= end_base) );

tic;
for i = 1:length(unique_coherence)
    for k = 1:2:3
        for g = 1:2
            %setup cols, y = mean psth, 
            %x: col1 = time, col2 = cuedir type, col3 = coherence*time(=slope), col4 = motion direc(=correct_vs_err),
            y_popramp = [y_popramp; mean_ungrouped_prestim_psth{i,1,k,g}(rampinds)'];
            y_poprampnorm = [y_poprampnorm; normmean_ungrouped_prestim_psth{i,1,k,g}(rampinds)'];
            x_popramp = [x_popramp; [0:end_ramp-begin_ramp]' [k.*ones(size(rampinds'))] [unique_coherence(i).*ones(size(rampinds'))].*[0:end_ramp-begin_ramp]' [unique_direction(g).*ones(size(rampinds'))] [unique_coherence(i).*ones(size(rampinds'))] ];
            y_popbase = [y_popbase; mean_ungrouped_prestim_psth{i,1,k,g}(baseinds)'];
            y_popbasenorm = [y_popbasenorm; normmean_ungrouped_prestim_psth{i,1,k,g}(baseinds)'];
            x_popbase = [x_popbase; [0:end_base-begin_base]' [k.*ones(size(baseinds'))] [unique_coherence(i).*ones(size(baseinds'))].*[0:end_base-begin_base]' [unique_direction(g).*ones(size(baseinds'))] [unique_coherence(i).*ones(size(baseinds'))] ];
            temp = gooddata;
            temp(excluded_cells_ungrouped_data{i,1,k,g}) = 0;
            temp = find(temp); %yields a list of indices for cells that are both in gooddata and have data for this condition
            for m = 1:length(temp)
                %setup cols, y = mean psth,
                %x: col1 = time, col2 = cuedir type, col3 = coherence*time(=slope), 
                %   col4 = motion direc(=correct_vs_err), col5 = cell#, col6 = zero_pct_xing (vs neutral for that cell)
                y_cellramp = [y_cellramp; cum_ungrouped_prestim_psth{i,1,k,g}(temp(m),rampinds)'];
                x_cellramp = [x_cellramp; [0:end_ramp-begin_ramp]' [k.*ones(size(rampinds'))] [unique_coherence(i).*ones(size(rampinds'))].*[0:end_ramp-begin_ramp]' [unique_direction(g).*ones(size(rampinds'))] [temp(m).*ones(size(rampinds'))] [(cum_zero_pct_xing(temp(m),k)-cum_zero_pct_xing(temp(m),2)).*ones(size(rampinds'))] [unique_coherence(i).*ones(size(rampinds'))]];
                y_cellbase = [y_cellbase; cum_ungrouped_prestim_psth{i,1,k,g}(temp(m),baseinds)'];
                x_cellbase = [x_cellbase; [0:end_base-begin_base]' [k.*ones(size(baseinds'))] [unique_coherence(i).*ones(size(baseinds'))].*[0:end_base-begin_base]' [unique_direction(g).*ones(size(baseinds'))] [temp(m).*ones(size(baseinds'))] [(cum_zero_pct_xing(temp(m),k)-cum_zero_pct_xing(temp(m),2)).*ones(size(baseinds'))] [unique_coherence(i).*ones(size(baseinds'))]];
                if ~isnan(cum_ungrouped_prestim_psth{i,1,k,g}(temp(m),1)) && ~isnan(cum_ungrouped_prestim_psth{i,1,2,g}(temp(m),1))
                    slope = regress(cum_ungrouped_prestim_psth{i,1,k,g}(temp(m),rampinds)',[ones(length(rampinds),1) (1:length(rampinds))']);
                    slope_neu = regress(cum_ungrouped_prestim_psth{i,1,2,g}(temp(m),rampinds)',[ones(length(rampinds),1) (1:length(rampinds))']);
                else
                    slope = [NaN NaN]'; slope_neu = [NaN NaN]'; 
                end
                cellrampavg = [cellrampavg; (slope-slope_neu)' k unique_coherence(i) unique_direction(g) (cum_zero_pct_xing(temp(m),k)-cum_zero_pct_xing(temp(m),2)), temp(m)];
                temp2 = mean(cum_ungrouped_prestim_psth{i,1,k,g}(temp(m),baseinds)) - mean(cum_ungrouped_prestim_psth{i,1,2,g}(temp(m),baseinds));
                cellbaseavg = [cellbaseavg; temp2 k unique_coherence(i) unique_direction(g) (cum_zero_pct_xing(temp(m),k)-cum_zero_pct_xing(temp(m),2)) temp(m)];
            end
        end
    end
end
toc
disp('x: col1 = time, col2 = cuedir type, col3 = coherence*time(=slope), col4 = motion direc(=correct_vs_err)')
stepwisefit(x_popramp,y_popramp); %note that stepwise fit does NOT require a column of 1s for offset.          
stepwisefit(x_popbase,y_popbase);
disp('x: col1 = time, col2 = cuedir type, col3 = coherence*time(=slope),')
disp('   col4 = motion direc(=correct_vs_err), col5 = cell#, col6 = zero_pct_xing (vs neutral for that cell)')
stepwisefit(x_cellramp,y_cellramp); %note that stepwise fit does NOT require a column of 1s for offset.          
stepwisefit(x_cellbase,y_cellbase);


%% I'm not sure but i believe the following code is to organize the data
%into a format where it can be saved as a text file and easily imported
%into statistica for further analyses.
fullbasematrix = []; fullbasematrix_labels = [];
rowcounter = 1; colcounter = 1;
motdirlabels = {'mP','mN'}; cuedirlabels = {'cNul','cNeu','cPrf'};
for i = 1:length(unique_coherence)
    for k = 1:3
        for g = 1:2
            fullbasematrix_labels = strcat(fullbasematrix_labels, sprintf(' %d_%s_%s',unique_coherence(i),cuedirlabels{k},motdirlabels{g}));
        end
    end
end
fullbasematrix_labels = sprintf('%s\r\n',fullbasematrix_labels);
% fid = fopen('fullbasematrix.txt','a');
% fprintf(fid,fullbasematrix_labels);
% for m = 1:length(gooddata)
%     if (gooddata(m))
%         colcounter = 1;
%         for i = 1:length(unique_coherence)
%             for k = 1:3 %cuedir
%                 for g = 1:2 %motion dir
%                     if ~isempty(find(excluded_cells_ungrouped_data{i,1,k,g}==m))
%                         fullbasematrix(end+1,:) = [unique_coherence(i) k g -9999];
%                         slopematrix(end+1,:) = [unique_coherence(i) k g -9999];
%                     else
%                         fullbasematrix(end+1,:) = [unique_coherence(i) k g mean(cum_ungrouped_prestim_psth{i,1,k,g}(m,baseinds))];
%                         [slope intercept] = regress(cum_ungrouped_prestim_psth{i,1,k,g}(temp(m),rampinds)',(1:length(rampinds))');
%                         slopematrix(end+1,:) = [unique_coherence(i) k g slope];
%                     end
%                     colcounter = colcounter+1;
%                 end
%             end
%         end
%         fprintf(fid,sprintf('%s\r\n',num2str(fullbasematrix(rowcounter,:))));
%         rowcounter = rowcounter+1;
%     end
% end
% fclose(fid);
% save('basematrix.txt','-ascii','basematrix')
% for m = 1:length(gooddata)
%     if (gooddata(m))
%         colcounter = 1;
%         for i = 1:length(unique_coherence)
%             for k = 1:3 %cuedir
%                 for g = 1:2 %motion dir
%                     if ~isempty(find(excluded_cells_ungrouped_data{i,1,k,g}==m))
%                         fullbasematrix(rowcounter,colcounter) = -9999; %statistica code for missing element
%                     else
%                         fullbasematrix(rowcounter,colcounter) = mean(cum_ungrouped_prestim_psth{i,1,k,g}(m,baseinds));
% %                         [slope intercept] = regress(cum_ungrouped_prestim_psth{i,1,k,g}(temp(m),rampinds)',(1:length(rampinds))');
%                     end
%                     colcounter = colcounter+1;
%                 end
%             end
%         end
%         fprintf(fid,sprintf('%s\r\n',num2str(fullbasematrix(rowcounter,:))));
%         rowcounter = rowcounter+1;
%     end
% end

%% test for changes in offset for motion-aligned PSTHs
%first subtract neutral cue response for each cell. 
%second, average across all conditions (5 coherences * 2 choices)
%NB: matlab's ttest ignores those elements in the array with value NaN (so don't need to worry about cells that don't have a trials in a particular condition)
t_on = 245; t_off = 445;
comp = [1 2];% 2 3; 1 2]
for c = 1:size(comp,1)
    counter = 0; clear meanpsthdiff grouppsthdiff psthdiff;
    psthdiff = cell(length(unique_coherence),2,2);
    meanpsthdiff = cell(length(unique_coherence),2);
    grouppsthdiff = cell(2);
    for i = 1:length(unique_coherence)
        for j = 1:2 %choice
            counter = counter + 1;
            for kind = 1:2 %only operate on T1- and T2-cued trials
                k = comp(c,kind);
                psthdiff{i,j,kind} = cum_prestim_psth{i,j,k}(gooddata,prestim+t_on+1:prestim+t_off) - cum_prestim_psth{i,j,2}(gooddata,prestim+t_on+1:prestim+t_off);
                meanpsthdiff{i,j}(:,kind) = mean(psthdiff{i,j,kind},2); %averages across time in psth into delta firing rates (with funny units)
                grouppsthdiff{kind}(:,(i-1)*2+j) = meanpsthdiff{i,j}(:,kind);
                %reorganize the average firing rates so that each matrix has a
                %list of the firing rate differences for the 10 combinations of coherence and choice.  This is for subsequent testing.
            end
            [h_offset(i,j) p_offset(i,j)] = ttest(meanpsthdiff{i,j}(:,1), meanpsthdiff{i,j}(:,2));
        end
    end
    [h_offset p_offset]
    [h_groupoffset p_groupoffset] = ttest(mean(grouppsthdiff{1},2), mean(grouppsthdiff{2},2))
    [nanmean(nanmean(grouppsthdiff{1},2))-nanmean(nanmean(grouppsthdiff{2},2))]
end

%% timecourse of firing rate diff (relative to neutral), 
win_step = 25; win_halfwidth = 25;
win_ctr = [0:win_step:700];
clear meanratediff cohmeanratediff fullmeanratediff cellmeanratediff cellstdratediff
meanratediff = zeros(length(win_ctr),3,5,2,sum(gooddata));
fullmeanratediff = zeros(length(win_ctr),3,5,sum(gooddata));
cellmeanratediff = zeros(length(win_ctr),3,5);
cellstdratediff = zeros(length(win_ctr),3,5);
nonzero_meanratediff = zeros(length(win_ctr),3);
nonzero_stdratediff = zeros(length(win_ctr),3);
allmeanratediff = zeros(length(win_ctr),3);
allstdratediff = zeros(length(win_ctr),3);
for t = 1:length(win_ctr)
    for k = 1:3
        for i = 1:5
            for j = 1:2
                t_win = prestim+win_ctr(t)-win_halfwidth : prestim+win_ctr(t)+win_halfwidth;
                meanratediff(t,k,i,j,:) = nanmean(cum_prestim_psth{i,j,k}(gooddata,t_win)-cum_prestim_psth{i,j,2}(gooddata,t_win),2);
            end
            fullmeanratediff(t,k,i,:) = nanmean(meanratediff(t,k,i,:,:),4); %avg across choice
            cellmeanratediff(t,k,i) = nanmean(fullmeanratediff(t,k,i,:)); %avg across cells
            cellstdratediff(t,k,i) = nanstd(fullmeanratediff(t,k,i,:))./sqrt(sum(gooddata)); %std err across cells
        end
        nonzero_meanratediff(t,k) = nanmean(cellmeanratediff(t,k,2:5)); %avg across coher and cell and choice
        nonzero_stdratediff(t,k) = squeeze(nanstd(nanmean(fullmeanratediff(t,k,2:5,:),3)))./sqrt(sum(gooddata)); %cell-stderr of avg of coher and chioce
        allmeanratediff(t,k) = nanmean(cellmeanratediff(t,k,:)); %avg across coher and cell and choice
        allstdratediff(t,k) = squeeze(nanstd(nanmean(fullmeanratediff(t,k,:,:),3)))./sqrt(sum(gooddata)); %cell-stderr of avg of coher and chioce
    end
end
h = figure; hold on
for k = 1:3
    for i = 1:5
        subplot(6,1,i); hold on;
%         plot(win_ctr,cellmeanratediff(:,k,i),linesymbols{k});
        errorbar(win_ctr,cellmeanratediff(:,k,i),cellstdratediff(:,k,i),linesymbols{k});
        if k==1, ylabel(sprintf('Coh=%d%%',unique_coherence(i))); end
    end
    subplot(616); hold on;
%     plot(win_ctr,allmeanratediff(:,k),linetypes{k});
    errorbar(win_ctr,allmeanratediff(:,k),allstdratediff(:,k),linetypes{k});
    if k==1
        ylabel('All Coher'); xlabel('Time Post Motion Onset (ms)'); 
    end
end
temp = [win_ctr' reshape([nonzero_meanratediff; nonzero_stdratediff], length(win_ctr),6)];
% save('Z:\LabTools\Matlab\TEMPO_Analysis\ProtocolSpecific\CuedDirectionDiscrim\group_FRdiff_timecourse_shortdurdelay.txt', '-ascii', 'temp')

%% plot mean firing rate across a defined window, as a function of coherence, relative to neutral cue firing rate.
h = figure; hold on
t_on = startperiod; t_off = endperiod;
clear psthdiff meanratediff cohmeanratediff fullmeanratediff cellmeanratediff cellstdratediff 
for k = 1:3
    for i = 1:5
        for j=1:2
            psthdiff{i,j,k} = cum_ungrouped_prestim_psth{i,j,k,1}(gooddata,prestim+t_on+1:prestim+t_off) - cum_ungrouped_prestim_psth{i,j,2,1}(gooddata,prestim+t_on+1:prestim+t_off);
            meanratediff(k,i,j,:) = nanmean(psthdiff{i,j,k}, 2); %avg across time 
        end
        fullmeanratediff(k,i,:) = nanmean(meanratediff(k,i,:,:),3); %avg across choice
        cellmeanratediff(k,i) = nanmean(fullmeanratediff(k,i,:)); %avg across cells
        cellstdratediff(k,i) = nanstd(fullmeanratediff(k,i,:)); %std dev across cells
    end
    cohmeanratediff(k,:) = nanmean(fullmeanratediff(k,:,:),2); %avg across coher, separated by cell and cuedir
    errorbar(unique_coherence, cellmeanratediff(k,:), cellstdratediff(k,:)./sqrt(sum(gooddata)),linesymbols{k});
    for j = 1:2
        for i = 1:5
            choicemeanratediff(k,i,j) = nanmean(squeeze(meanratediff(k,i,j,:)));
            choicestdratediff(k,i,j) = nanstd(squeeze(meanratediff(k,i,j,:)));
        end
%         errorbar(unique_coherence, choicemeanratediff(k,:,j), choicestdratediff(k,:,j)./sqrt(sum(gooddata)),linesymbolsT{j,k});
    end
end
xlabel('Coherence');
ylabel('Delta Firing Rate');
title(sprintf('time on = %d, time off = %d',t_on, t_off));
temp = [unique_coherence' cellmeanratediff(1,:)' cellstdratediff(1,:)'./sqrt(sum(gooddata)) ...
    cellmeanratediff(2,:)' cellstdratediff(2,:)'./sqrt(sum(gooddata)) cellmeanratediff(3,:)' cellstdratediff(3,:)'./sqrt(sum(gooddata))];
% save(sprintf('Z:\LabTools\Matlab\TEMPO_Analysis\ProtocolSpecific\CuedDirectionDiscrim\group_delFRxCoh_shortdurdelay%d%d.txt',t_on,t_off), '-ascii', 'temp')
cohmeanratediff(1,:)'-cohmeanratediff(2,:)';%null-neu for each cell - unsemicolon to view and copy
cohmeanratediff(3,:)'-cohmeanratediff(2,:)';%pref-neu for each cell
cohmeanratediff(3,:)'-cohmeanratediff(1,:)';%pref-null for each cell
        
%% test for changes in slope for motion-aligned PSTHs
t_on = startperiod; t_off = endperiod-50;
for i = 1:length(unique_coherence)
    counter = 0;
    for m = 1:length(cum_file)
        if gooddata(m)
            counter = counter+1;
            for k = 1:3
                temp = cum_ungrouped_prestim_psth{i,1,k,1}(m,prestim+t_on+1:prestim+t_off)';
                if isnan(temp(1))
                    cue_slope{i,k}(counter) = NaN;                    
                else
                    b_temp = regress(temp, [ones(size(temp)) [1:t_off-t_on]']);
                    cue_slope{i,k}(counter) = b_temp(2);
                end
            end
        end
    end
    diff_slope{i}(:,1) = cue_slope{i,3}-cue_slope{i,2}; %first column is pref-neutral slope for each neuron; cells are coherences
    diff_slope{i}(:,2) = cue_slope{i,1}-cue_slope{i,2}; %second column is null-neutral slope
    [h_slope(i) p_slope(i)] = ttest(diff_slope{i}(:,1), diff_slope{i}(:,2));
    groupdiff_slope{1}(i,:) = diff_slope{i}(:,1)'; %within each cell, each row has different neurons' slope for one coherence
    groupdiff_slope{2}(i,:) = diff_slope{i}(:,2)';
    for k = 1:3
        meanslope(i,k) = nanmean(cue_slope{i,k}); 
        seslope(i,k) = nanstd(cue_slope{i,k})./sqrt(sum(~isnan(cue_slope{i,k})));
    end
end
[h_slope' p_slope']
meandiff_slope = [nanmean(groupdiff_slope{1},1)' nanmean(groupdiff_slope{2},1)'];
[h_groupslope p_groupslope] = ttest(meandiff_slope(:,1), meandiff_slope(:,2))
mean(meandiff_slope,1)
plotmeandiffslope{1} = nanmean(groupdiff_slope{1},2); plotsediffslope{1} = nanstd(groupdiff_slope{1}')./sqrt(sum(~isnan(sum(groupdiff_slope{1}+groupdiff_slope{1}))));
plotmeandiffslope{2} = nanmean(groupdiff_slope{2},2); plotsediffslope{2} = nanstd(groupdiff_slope{2}')./sqrt(sum(~isnan(sum(groupdiff_slope{1}+groupdiff_slope{2}))));
figure; hold on;    
errorbar(unique_coherence, plotmeandiffslope{1}, plotsediffslope{1},'g>-');
errorbar(unique_coherence, plotmeandiffslope{2}, plotsediffslope{2},'b*-');
plot(unique_coherence, [0 0 0 0 0], 'ro-');
xlabel('Coherence (%)'); ylabel('Slope (Hz/ms)')
temp = [unique_coherence' 1000.*plotmeandiffslope{2} 1000.*plotsediffslope{2}' zeros(5,1) zeros(5,1) ...
    1000.*plotmeandiffslope{1} 1000.*plotsediffslope{1}' 1000.*meanslope 1000.*seslope];
% save(sprintf('Z:\\LabTools\\Matlab\\TEMPO_Analysis\\ProtocolSpecific\\CuedDirectionDiscrim\\group_slopexCoh_shortdurdelay%d%d.txt',t_on,t_off), '-ascii', 'temp')
%do 1-way repeated measures anova on neutral cue slope
y = [cue_slope{1,2} cue_slope{2,2} cue_slope{3,2} cue_slope{4,2} cue_slope{5,2}]';
cells = repmat((1:n)',5,1);
f1 = repmat((1:5)',n,1);
rmaov1([y f1 cells],0.05)
%do 2-way repeated measures anova on neutral cue slope with cue dir and
%coherence as factors
y = []; cells = []; f1 = []; f2 = [];
for k = 1:3
    y = [y cue_slope{1,2} cue_slope{2,2} cue_slope{3,2} cue_slope{4,2} cue_slope{5,2}];
    cells = [cells; repmat((1:n)',5,1)];
    f1 = [f1; repmat((1:5)',n,1)]; %coherence
    f2 = [f2; k.*ones(5*n,1)]; %cue dir
end
rmaov2([y' f1 f2 cells])
% f1name = 'Coher';f2name = 'Cuedir';
% rm_anova2(y,cells,f1,f2,{f1name,f2name})
%% view the uptick on a cell by cell basis
figure; hold on;
counter = 0; map = colormap;
for m = 1:length(cum_file)
    if gooddata(m)
        counter = counter+1;
        temp = mean([cum_ungrouped_prestim_psth{1,1,1,1}(m,501:701); cum_ungrouped_prestim_psth{2,1,1,1}(m,501:701);...
            cum_ungrouped_prestim_psth{3,1,1,1}(m,501:701); cum_ungrouped_prestim_psth{4,1,1,1}(m,501:701);...
            cum_ungrouped_prestim_psth{5,1,1,1}(m,501:701)]);
        plot(prestim_x(501:701), temp, 'Color', map(ceil(counter/sum(gooddata)*64),:));
    end
end
for k = 1:3
    temp = mean([mean_ungrouped_prestim_psth{1,1,k,1}(501:701); mean_ungrouped_prestim_psth{2,1,k,1}(501:701);...
        mean_ungrouped_prestim_psth{3,1,k,1}(501:701); mean_ungrouped_prestim_psth{4,1,k,1}(501:701);...
        mean_ungrouped_prestim_psth{5,1,k,1}(501:701)]);
    plot(prestim_x(501:701),temp,linetypes{k},'linewidth',3);
end

%% compute CMI (Cue Modulation Index) for Valid and Invalid trials for each cell and plot scatter
%CMI = (cue - neutral) / (cue + neutral), taking the mean response during the specified window
t_on = 220; t_off = 420;
clear CMI meanCMI p;
for i = 1:length(unique_coherence)
    for j = 1:2
        for k = 1:3
            tempmean{i,j}(k,:) = mean(cum_prestim_psth{i,j,k}(:,prestim+t_on+1:prestim+t_off),2);
        end
        CMI(i,j,1,:) = (tempmean{i,j}(1,gooddata) - tempmean{i,j}(2,gooddata)) ./ ...
                       (tempmean{i,j}(1,gooddata) + tempmean{i,j}(2,gooddata));
        CMI(i,j,2,:) = (tempmean{i,j}(3,gooddata) - tempmean{i,j}(2,gooddata)) ./ ...
                       (tempmean{i,j}(3,gooddata) + tempmean{i,j}(2,gooddata));
    end
end
meanCMI(1,:) = nanmean(nanmean(CMI(:,:,1,:),1),2); %avg across coherence and choice
meanCMI(2,:) = nanmean(nanmean(CMI(:,:,2,:),1),2);            
[h,p(1)] = ttest(meanCMI(1,:),meanCMI(2,:));
[mean(meanCMI,2) std(meanCMI,0,2)]
figure;
plot(meanCMI(2,:),meanCMI(1,:),'b*')
xlabel('T2-cue CMI'); ylabel('T1-cue CMI');
title(sprintf('CMI: %dms to %dms post-motion onset', t_on, t_off));
%xlim([-.2,.2]); ylim([-.2,.2]);
[h,p(2)] = ttest(meanCMI(2,:),0);
[h,p(3)] = ttest(meanCMI(1,:),0);
[rho,corrp] = corr(meanCMI');
[p corrp(2) rho(2)]
% keyboard;

                
%% plot estimates of latency
% figure;
% offset = 0;
% for i = 1:length(unique_coherence)
%     subplot(5,1,i); hold on
%     for k = 1:3
%         diff_meanpsth = mean_prestim_psth{i,1,k} - mean_prestim_psth{i,2,k};
%         t = 0; keep_looping = 1;
%         while keep_looping
%             preinds = find(prestim_x < t);
%             thresh = mean(diff_meanpsth(preinds)) + 3.*std(diff_meanpsth(preinds));
%             if ( (diff_meanpsth(prestim+t) > thresh) | (t == max(prestim_x)) );
%                 keep_looping = 0;
%                 mean_lat(i,k) = t + offset;
%             else
%                 t = t+1;
%             end
%         end        
%         plot(prestim_x,diff_meanpsth,linetypes{k}, 'LineWidth',1.2);
%     end
%     axis tight;
%     plot(repmat(mean_lat(i,1),2,1),ylim,linetypes{1},'LineWidth',1.2);
%     plot(repmat(mean_lat(i,2),2,1),ylim,linetypes{2},'LineWidth',1.2);
%     plot(repmat(mean_lat(i,3),2,1),ylim,linetypes{3},'LineWidth',1.2);
% end
% keyboard

%% plot choice-sorted neutral cue psths on same axes at different coherences
h = figure ; hold on;
set(h,'PaperPosition', [.2 .2 8 10.7], 'Position', [250 50 500 573], 'Name', 'PSTH by Coherence');
temp = prestim_x';
for k = 1:3
    legstr = 'legh(k) = legend(coh_h';
    subplot(3,1,k); hold on;
    for i = 1:length(unique_coherence)
%         coh_h(i) = plot(prestim_x, mean_ungrouped_prestim_psth{i,1,k,1},coh_lines{i});
        coh_h(i) = plot(prestim_x, mean_ungrouped_prestim_psth{i,1,k,1},coh_lines{i},'LineWidth',1.2);
        plot(prestim_x, mean_ungrouped_prestim_psth{i,2,k,2},coh_dashedlines{i},'LineWidth',1.2);
        legstr = strcat(legstr, sprintf(',''%d%%''',unique_coherence(i)));
    end
    axis tight;
    legstr = strcat(legstr, ',''Location'',''NorthWest'');');
    eval(legstr); set(legh(k),'box','off');
    title(sprintf('Correct PrefChoice Trials: %s',unique_cue_dir_type_names{k}));
    ylabel('Firing Rate (Hz)');
    if k==3 
        xlabel('Time about Motion Onset'); 
    end
end
temp = prestim_x';
for i = 1:5, temp = [temp mean_prestim_psth{i,1,2}']; end
for i = 1:5, temp = [temp mean_prestim_psth{i,2,2}']; end
% save('Z:\LabTools\Matlab\TEMPO_Analysis\ProtocolSpecific\CuedDirectionDiscrim\group_prestimPSTHxCoh_shortdurdelay.txt', '-ascii', 'temp')
temp = sacc_x';
for i = 1:5, temp = [temp mean_sacc_psth{i,1,2}']; end
for i = 1:5, temp = [temp mean_sacc_psth{i,2,2}']; end
% save('Z:\LabTools\Matlab\TEMPO_Analysis\ProtocolSpecific\CuedDirectionDiscrim\group_saccPSTHxCoh_shortdurdelay.txt', '-ascii', 'temp')
%% plot motion-sorted neutral cue psths on same axes at different coherences
h = figure ; hold on;
set(h,'PaperPosition', [.2 .2 8 10.7], 'Position', [250 50 500 573], 'Name', 'PSTH by Coherence');
coh_lines = {'r','m','b','c','g'};
coh_dashedlines = {'r--','m--','b--','c--','g--'};
coh_dottedlines = {'r:','m:','b:','c:','g:'};
temp = prestim_x';
for k = 1:3
    legstr = 'legh(k) = legend(coh_h';
    subplot(3,1,k); hold on;
    for i = 1:length(unique_coherence)
%         coh_h(i) = plot(prestim_x, mean_ungrouped_prestim_psth{i,1,k,1},coh_lines{i});
        coh_h(i) = plot(prestim_x, mean_motgrouped_prestim_psth{i,k,1},coh_lines{i},'LineWidth',1.2);
        plot(prestim_x, mean_motgrouped_prestim_psth{i,k,2},coh_dashedlines{i},'LineWidth',1.2);
        legstr = strcat(legstr, sprintf(',''%d%%''',unique_coherence(i)));
    end
    axis tight;
    legstr = strcat(legstr, ',''Location'',''NorthWest'');');
    eval(legstr); set(legh(k),'box','off');
    title(sprintf('Motion-sorted Trials: %s',unique_cue_dir_type_names{k}));
    ylabel('Firing Rate (Hz)');
    if k==3 
        xlabel('Time about Motion Onset'); 
    end
end
temp = prestim_x';
for i = 1:5, temp = [temp mean_motgrouped_prestim_psth{i,2,1}']; end
for i = 1:5, temp = [temp mean_motgrouped_prestim_psth{i,2,2}']; end
% save('Z:\LabTools\Matlab\TEMPO_Analysis\ProtocolSpecific\CuedDirectionDiscrim\group_motsort_prestimPSTHxCoh_shortdurdelay.txt', '-ascii', 'temp')


%% %now perform paired t-tests comparing prefdir cues and null cues
tic
comp = [1 3];
comp_co = [1 3]; %FIXED... no neutral in cue only
comp_color = {'b','g'}; 
comp_color_co = {'b','g'};
binwidth = 21; 
tstep = 1;
for v=1:tstep:length(postcue_x)
    [h_postcue_combined(v) p_postcue_combined(v)] = ...
        ttest(cum_postcue_combined_psth{comp(1)}(gooddata,v), cum_postcue_combined_psth{comp(2)}(gooddata,v),0.05,'both');
%     [h_norm_postcue_combined(v) p_norm_postcue_combined(v)] = ...
%         ttest(norm_postcue_combined_psth{comp(1)}(gooddata,v), norm_postcue_combined_psth{comp(2)}(gooddata,v),0.05,'both');
    for j = 1:2
        [h_postcue{j}(v) p_postcue{j}(v)] = ...
            ttest(cum_postcue_psth{j,comp(1)}(gooddata,v), cum_postcue_psth{j,comp(2)}(gooddata,v),0.05,'both');
%         [h_norm_postcue{j}(v) p_norm_postcue{j}(v)] = ...
%             ttest(norm_postcue_psth{j,comp(1)}(gooddata,v), norm_postcue_psth{j,comp(2)}(gooddata,v),0.05,'both');
    end
end
for v=1:length(prestim_x)
    for j = 1:2
        for i = 1:length(unique_coherence)
                [h_prestim{i,j}(v) p_prestim{i,j}(v)] = ...
                ttest(cum_prestim_psth{i,j,comp(1)}(gooddata,v), cum_prestim_psth{i,j,comp(2)}(gooddata,v),0.05,'both');
%             [h_norm_prestim{i,j}(v) p_norm_prestim{i,j}(v)] = ...
%                 ttest(norm_prestim_psth{i,j,comp(1)}(gooddata,v), norm_prestim_psth{i,j,comp(2)}(gooddata,v),0.05,'both');
                [h_lort_prestim{i,j}(v) p_lort_prestim{i,j}(v)] = ...
                ttest(cum_lort_prestim_psth{i,j,comp(1)}(gooddata,v), cum_lort_prestim_psth{i,j,comp(2)}(gooddata,v),0.05,'both');
                [h_hirt_prestim{i,j}(v) p_hirt_prestim{i,j}(v)] = ...
                ttest(cum_hirt_prestim_psth{i,j,comp(1)}(gooddata,v), cum_hirt_prestim_psth{i,j,comp(2)}(gooddata,v),0.05,'both');
                [h_loisi_prestim{i,j}(v) p_loisi_prestim{i,j}(v)] = ...
                ttest(cum_loisi_prestim_psth{i,j,comp(1)}(gooddata,v), cum_loisi_prestim_psth{i,j,comp(2)}(gooddata,v),0.05,'both');
                [h_hiisi_prestim{i,j}(v) p_hiisi_prestim{i,j}(v)] = ...
                ttest(cum_hiisi_prestim_psth{i,j,comp(1)}(gooddata,v), cum_hiisi_prestim_psth{i,j,comp(2)}(gooddata,v),0.05,'both');
            for g = 1:2
                [h_prestim_ungrouped{i,j,g}(v) p_prestim_ungrouped{i,j,g}(v)] = ...
                ttest(cum_ungrouped_prestim_psth{i,j,comp(1),g}(gooddata,v), cum_ungrouped_prestim_psth{i,j,comp(2),g}(gooddata,v),0.05,'both');
            end
        end
    end
end
for v=1:length(sacc_x)
    for j = 1:2
        [h_co_sacc{j}(v) p_co_sacc{j}(v)] = ...
            ttest(cum_co_sacc_psth{j,comp_co(1)}(gooddata_co,v), cum_co_sacc_psth{j,comp_co(2)}(gooddata_co,v),0.05,'both');
%         [h_norm_co_sacc{j}(v) p_norm_co_sacc{j}(v)] = ...
%             ttest(norm_co_sacc_psth{j,comp(1)}(gooddata_co,v), norm_co_sacc_psth{j,comp(2)}(gooddata_co,v),0.05,'both');
        for i = 1:length(unique_coherence)
            [h_sacc{i,j}(v) p_sacc{i,j}(v)] = ...
                ttest(cum_sacc_psth{i,j,comp(1)}(gooddata,v), cum_sacc_psth{i,j,comp(2)}(gooddata,v),0.05,'both');
%             [h_norm_sacc{i,j}(v) p_norm_sacc{i,j}(v)] = ...
%                 ttest(norm_sacc_psth{i,j,comp(1)}(gooddata,v), norm_sacc_psth{i,j,comp(2)}(gooddata,v),0.05,'both');
        end
    end
end
for v=1:tstep:length(postcue_co_x)
    for j = 1:2
        [h_co_postcue{j}(v) p_co_postcue{j}(v)] = ...
            ttest(cum_co_postcue_psth{j,comp_co(1)}(gooddata_co,v), cum_co_postcue_psth{j,comp_co(2)}(gooddata_co,v),0.05,'both');
%         [h_norm_co_postcue{j}(v) p_norm_co_postcue{j}(v)] = ...
%             ttest(norm_co_postcue_psth{j,comp(1)}(gooddata_co,v), norm_co_postcue_psth{j,comp(2)}(gooddata_co,v),0.05,'both');
    end
end
for v=1:length(prestim_co_x)
    for j = 1:2
        [h_co_prestim{j}(v) p_co_prestim{j}(v)] = ...
            ttest(cum_co_prestim_psth{j,comp_co(1)}(gooddata_co,v), cum_co_prestim_psth{j,comp_co(2)}(gooddata_co,v),0.05,'both');
%         [h_norm_co_prestim{j}(v) p_norm_co_prestim{j}(v)] = ...
%             ttest(norm_co_prestim_psth{j,comp(1)}(gooddata_co,v), norm_co_prestim_psth{j,comp(2)}(gooddata_co,v),0.05,'both');
    end
end
toc

%% now do permutation test, shuffling across assignment of two cuedirs in comp
% nboot = 1000; alpha = 0.05;
% CIlobnd = floor(nboot*(alpha/2)); CIupbnd = ceil(nboot*(1-alpha/2)); %demarcate central 95%
% for b = 1:nboot
%     goodcells = find(gooddata);
%     for i = 1:length(goodcells)
%         temp = ceil(2*rand);
%         boot_cum_psth{1}(i,:) = cum_prestim_psth{comp(temp)}(goodcells(i),:);
%         boot_cum_psth{2}(i,:) = cum_prestim_psth{comp(3-temp)}(goodcells(i),:);
%     end
%     boot_psth{1}(b,:) = nanmean(boot_cum_psth{1},1);
%     boot_psth{2}(b,:) = nanmean(boot_cum_psth{2},1);
% end
% for v = 1:length(prestim_x)
%     temp = sort(boot_psth{1}(:,v));
%     prestim_psth_CI{1}(:,v) = [temp(CIlobnd); temp(CIupbnd)];
%     temp = sort(boot_psth{2}(:,v));
%     prestim_psth_CI{2}(:,v) = [temp(CIlobnd); temp(CIupbnd)];
% end


%% anovas to compare activity across timecourses... 
% p_anova_postcue_combined = anova2([cum_postcue_combined_psth{1}(gooddata,:); cum_postcue_combined_psth{2}(gooddata,:); cum_postcue_combined_psth{3}(gooddata,:)], n,'off');
% for j = 1:2
%     p_anova_postcue(j,:) = anova2([cum_postcue_psth{j,1}(gooddata,:); cum_postcue_psth{j,2}(gooddata,:); cum_postcue_psth{j,3}(gooddata,:)], n,'off');
% %     p_anova_co_postcue(j,:) = anova2([cum_co_postcue_psth{j,1}(gooddata_co,:); cum_co_postcue_psth{j,3}(gooddata_co,:)], n_co,'off');
% %     p_anova_co_prestim(j,:) = anova2([cum_co_prestim_psth{j,1}(gooddata_co,:); cum_co_prestim_psth{j,3}(gooddata_co,:)], n_co,'off');
% %     p_anova_co_sacc(j,:) = anova2([cum_co_sacc_psth{j,1}(gooddata_co,:); cum_co_sacc_psth{j,3}(gooddata_co,:)], n_co,'off');
%     for i = 1:length(coher)
%         p_anova_prestim(i,j,:) = anova2([cum_prestim_psth{i,j,1}(gooddata,:); cum_prestim_psth{i,j,2}(gooddata,:); cum_prestim_psth{i,j,3}(gooddata,:)], n,'off');
%         p_anova_sacc(i,j,:) = anova2([cum_sacc_psth{i,j,1}(gooddata,:); cum_sacc_psth{i,j,2}(gooddata,:); cum_sacc_psth{i,j,3}(gooddata,:)], n,'off');
%     end
% end
% %repeated measure anovas (type help rm_anova2 for explanation of variables)
% y = [reshape(cum_postcue_combined_psth{1}(gooddata,:)',prod(size(cum_postcue_combined_psth{1}(gooddata,:))),1); reshape(cum_postcue_combined_psth{2}(gooddata,:)',prod(size(cum_postcue_combined_psth{2}(gooddata,:))),1); reshape(cum_postcue_combined_psth{3}(gooddata,:)',prod(size(cum_postcue_combined_psth{3}(gooddata,:))),1)];
% s = reshape(repmat([1:n], size(cum_postcue_combined_psth{1},2), 3),size(y));
% f1 = repmat([1:size(cum_postcue_combined_psth{1},2)]', length(y)/size(cum_postcue_combined_psth{1},2),1);
% f2 = reshape(repmat([1:3], prod(size(cum_postcue_combined_psth{1}(gooddata,:))),1),size(y));
% f1name = 'Time';  f2name = 'CueDirType';
% % p_rm_anova_postcue_combined_full = rm_anova2(y,s,f1,f2,{f1name, f2name});
% % p_rm_anova_postcue_combined = [p_rm_anova_postcue_combined{3,6} p_rm_anova_postcue_combined_full{4,6}]; %[p_cuedirtype p_time_X_cuedirtype]
% for j = 1:2
%     y = [reshape(cum_postcue_psth{j,1}(gooddata,:)',prod(size(cum_postcue_psth{j,1}(gooddata,:))),1); reshape(cum_postcue_psth{j,2}(gooddata,:)',prod(size(cum_postcue_psth{j,2}(gooddata,:))),1); reshape(cum_postcue_psth{j,3}(gooddata,:)',prod(size(cum_postcue_psth{j,3}(gooddata,:))),1)];
%     s = reshape(repmat([1:n], size(cum_postcue_psth{j,1},2), 3),size(y));
%     f1 = repmat([1:size(cum_postcue_psth{j,1},2)]', length(y)/size(cum_postcue_psth{j,1},2), 1);
%     f2 = reshape(repmat([1:3], prod(size(cum_postcue_psth{j,1}(gooddata,:))),1),size(y));
% %     p_rm_anova_postcue_full{j} = rm_anova2(y,s,f1,f2,{f1name, f2name});
% %     p_rm_anova_postcue(j,:) = [p_rm_anova_postcue_full{j}{3,6} p_rm_anova_postcue_full{j}{4,6}]
%     for i = 1:length(coher)
%         y = [reshape(cum_prestim_psth{i,j,1}(gooddata,:)',prod(size(cum_prestim_psth{i,j,1}(gooddata,:))),1); reshape(cum_prestim_psth{i,j,2}(gooddata,:)',prod(size(cum_prestim_psth{i,j,2}(gooddata,:))),1); reshape(cum_prestim_psth{i,j,3}(gooddata,:)',prod(size(cum_prestim_psth{i,j,3}(gooddata,:))),1)];
%         s = reshape(repmat([1:n], size(cum_prestim_psth{i,j,1},2), 3),size(y));
%         f1 = repmat([1:size(cum_prestim_psth{i,j,1},2)]', length(y)/size(cum_prestim_psth{i,j,1},2), 1);
%         f2 = reshape(repmat([1:3], prod(size(cum_prestim_psth{i,j,1}(gooddata,:))),1),size(y));
% %         p_rm_anova_prestim_full{i,j} = rm_anova2(y,s,f1,f2,{f1name, f2name});
% %         p_rm_anova_prestim(i,j,:) = [p_rm_anova_prestim_full{i,j}{3,6} p_rm_anova_prestim_full{i,j}{4,6}]
%         
%         y = [reshape(cum_sacc_psth{i,j,1}(gooddata,:)',prod(size(cum_sacc_psth{i,j,1}(gooddata,:))),1); reshape(cum_sacc_psth{i,j,2}(gooddata,:)',prod(size(cum_sacc_psth{i,j,2}(gooddata,:))),1); reshape(cum_sacc_psth{i,j,3}(gooddata,:)',prod(size(cum_sacc_psth{i,j,3}(gooddata,:))),1)];
%         s = reshape(repmat([1:n], size(cum_sacc_psth{i,j,1},2), 3),size(y));
%         f1 = repmat([1:size(cum_sacc_psth{i,j,1},2)]', length(y)/size(cum_sacc_psth{i,j,1},2), 1);
%         f2 = reshape(repmat([1:3], prod(size(cum_sacc_psth{i,j,1}(gooddata,:))),1),size(y));
% %         p_rm_anova_sacc_full{i,j} = rm_anova2(y,s,f1,f2,{f1name, f2name});
% %         p_rm_anova_sacc(i,j,:) = [p_rm_anova_sacc_full{i,j}{3,6} p_rm_anova_sacc_full{i,j}{4,6}]
%     end
% 	y = [reshape(cum_co_postcue_psth{j,1}(gooddata_co,:)',prod(size(cum_co_postcue_psth{j,1}(gooddata_co,:))),1); reshape(cum_co_postcue_psth{j,3}(gooddata_co,:)',prod(size(cum_co_postcue_psth{j,3}(gooddata_co,:))),1)];
%     s = reshape(repmat([1:n_co], size(cum_co_postcue_psth{1},2), 2),size(y));
%     f1 = repmat([1:size(cum_co_postcue_psth{1},2)]', length(y)/size(cum_postcue_psth{1},2),1);
%     f2 = reshape(repmat([1 2], prod(size(cum_co_postcue_psth{1}(gooddata_co,:))),1),size(y));
% %     p_rm_anova_co_postcue_full{j} = rm_anova2(y,s,f1,f2,{f1name, f2name});
% %     p_rm_anova_co_postcue(j,:) = [p_rm_anova_co_postcue_full{j}{3,6} p_rm_anova_co_postcue_full{j}{4,6}]
%     
% 	y = [reshape(cum_co_prestim_psth{j,1}(gooddata_co,:)',prod(size(cum_co_prestim_psth{j,1}(gooddata_co,:))),1); reshape(cum_co_prestim_psth{j,3}(gooddata_co,:)',prod(size(cum_co_prestim_psth{j,3}(gooddata_co,:))),1)];
%     s = reshape(repmat([1:n_co], size(cum_co_prestim_psth{1},2), 2),size(y));
%     f1 = repmat([1:size(cum_co_prestim_psth{1},2)]', length(y)/size(cum_prestim_psth{1},2),1);
%     f2 = reshape(repmat([1 2], prod(size(cum_co_prestim_psth{1}(gooddata_co,:))),1),size(y));
% %     p_rm_anova_co_prestim_full{j} = rm_anova2(y,s,f1,f2,{f1name, f2name});
% %     p_rm_anova_co_prestim(j,:) = [p_rm_anova_co_prestim_full{j}{3,6} p_rm_anova_co_prestim_full{j}{4,6}]
%     
% 	y = [reshape(cum_co_sacc_psth{j,1}(gooddata_co,:)',prod(size(cum_co_sacc_psth{j,1}(gooddata_co,:))),1); reshape(cum_co_sacc_psth{j,3}(gooddata_co,:)',prod(size(cum_co_sacc_psth{j,3}(gooddata_co,:))),1)];
%     s = reshape(repmat([1:n_co], size(cum_co_sacc_psth{1},2), 2),size(y));
%     f1 = repmat([1:size(cum_co_sacc_psth{1},2)]', length(y)/size(cum_sacc_psth{1},2),1);
%     f2 = reshape(repmat([1 2], prod(size(cum_co_sacc_psth{1}(gooddata_co,:))),1),size(y));
% %     p_rm_anova_co_sacc_full{j} = rm_anova2(y,s,f1,f2,{f1name, f2name});
% %     p_rm_anova_co_sacc(j,:) = [p_rm_anova_co_sacc_full{j}{3,6} p_rm_anova_co_sacc_full{j}{4,6}]
% end    

% keyboard


%% now plot these averaged raw data
plot_ttests = 1;

%% first this figure contains motion psths for cued trials
handl(1) = figure; 
set(handl(1),'PaperPosition', [.2 .2 8 10.7], 'Position', [250 50 500 573], 'Name', 'Group data: Raw Motion PSTH');
for i = 1:length(unique_coherence)
    for j = 1:2
        subplot(length(unique_coherence),2, 2*(i-1)+j); hold on
        for k = 1:length(unique_cue_dir_type)
            plot(prestim_x, mean_prestim_psth{i,j,k}, linetypes{k}, 'LineWidth',1.2);
        end
        axis('tight'); upyl(j) = max(ylim); ylim([0 upyl(j)]);
        if j == 1
            ylabel(sprintf('Coh = %3.1f%%',unique_coherence(i)));
        end
        if i == 1
            if j==1
                title ('Pref Choices','Color','k');
            else
                title ('Null Choices','Color','k');
            end
        end
        if i == length(unique_coherence)
            xlabel ('Time about Motion onset');
        end
    end
    [maxyl, maxylind] = max(upyl);
    for j=1:2
        subplot(length(unique_coherence),2,2*(i-1)+j);
        if j~=maxylind
            ylim([0 maxyl]);
        end
        if plot_ttests
            if ~sum(isnan(h_prestim{i,j}))
                c1 = prestim_x(logical(h_prestim{i,j} & (mean_prestim_psth{i,j,comp(1)} > mean_prestim_psth{i,j,comp(2)})));
                c2 = prestim_x(logical(h_prestim{i,j} & (mean_prestim_psth{i,j,comp(1)} < mean_prestim_psth{i,j,comp(2)})));
                plot(c1, max(ylim).*ones(length(c1),1),sprintf('%s+',comp_color{1}),'MarkerSize',3);
                plot(c2, max(ylim).*ones(length(c2),1),sprintf('%s+',comp_color{2}),'MarkerSize',3);
            end
        end
        set(gca,'tickdir','out','xcolor','k','ycolor','k');
        fill([greybar_start greybar_end greybar_end greybar_start],[min(ylim) min(ylim) max(ylim) max(ylim)], [0.75 0.75 0.75], 'FaceAlpha',0.5,'LineStyle','none');
        fill([yellowbar_start yellowbar_end yellowbar_end yellowbar_start],[min(ylim) min(ylim) max(ylim) max(ylim)], [1 1 0], 'FaceAlpha',0.5,'LineStyle','none');
        line([endperiod endperiod],ylim,'Color','k','LineStyle',':');
        line([startperiod startperiod],ylim,'Color','k','LineStyle',':');
    end
end

%this formats labels and data for easy import into origin.  note that this
%does NOT import the t-test marker locations.
% temp = prestim_x';
% ilabels = {'0','2','4','8','16'};
% jlabels = {'T1','T2'};
% klabels = {'NLc','NUc','PRc'};
% templabels = 'Time ';
% for j = 1:2
%     for i = 1:5
%         for k = 1:3
%             temp = [temp mean_prestim_psth{i,j,k}'];
%             templabels = [templabels, ilabels{i}, jlabels{j}, klabels{k}, ' '];
%         end
%     end
% end
% fid = fopen('Z:\LabTools\Matlab\TEMPO_Analysis\ProtocolSpecific\CuedDirectionDiscrim\group_prestimPSTH_longdur.txt','W');
% fprintf(fid,templabels);
% fprintf(fid,'\r\n');
% for i = 1:size(temp,1)
%     fprintf(fid, '%9.5f\t', temp(i,:));
%     fprintf(fid, '\r\n');
% end
% fclose(fid)

%this is another version of the above, this time averaging across non-0 coherences and saving that
temp = prestim_x';
ilabels = {'0','n0'};
jlabels = {'T1','T2'};
klabels = {'NLc','NUc','PRc'};
templabels = 'Time ';
for j = 1:2
    for i = 1:2
        for k = 1:3
            if i==1 
                temp = [temp mean_prestim_psth{i,j,k}'];
            else %for non zero coherences, average across the non-zero psths and append that to temp
                crossmean_psth{j,k} = mean([mean_prestim_psth{2,j,k}; mean_prestim_psth{3,j,k}; ...
                      mean_prestim_psth{4,j,k}; mean_prestim_psth{5,j,k}], 1); 
                temp = [temp crossmean_psth{j,k}'];
            end
            templabels = [templabels, ilabels{i}, jlabels{j}, klabels{k}, ' '];
        end
    end
end
figure; hold on;
maxylim = 0;
for j = 1:2
    subplot(2,2,j); hold on; %zero-coh
    fill([greybar_start greybar_end greybar_end greybar_start],[0 0 50 50], [0.75 0.75 0.75], 'FaceAlpha',0.5,'LineStyle','none');
    fill([yellowbar_start yellowbar_end yellowbar_end yellowbar_start],[0 0 50 50], [1 1 0], 'FaceAlpha',0.5,'LineStyle','none');
    for k = 1:3, plot(prestim_x, mean_prestim_psth{1,j,k}, linetypes{k}, 'LineWidth',2); end
    if j==1, ylabel('Firing Rate (Hz)','Color','k'); end
    if max(ylim) > maxylim, maxylim = max(ylim); end
    subplot(2,2,2+j); hold on; %nonzero-coh
    fill([greybar_start greybar_end greybar_end greybar_start],[0 0 50 50], [0.75 0.75 0.75], 'FaceAlpha',0.5,'LineStyle','none');
    fill([yellowbar_start yellowbar_end yellowbar_end yellowbar_start],[0 0 50 50], [1 1 0], 'FaceAlpha',0.5,'LineStyle','none');
    for k = 1:3, plot(prestim_x, crossmean_psth{j,k}, linetypes{k}, 'LineWidth',2); end
    if j==1, ylabel('Firing Rate (Hz)','Color','k'); end
    xlabel('Time about Motion Onset (ms)','Color','k');
    if max(ylim) > maxylim, maxylim = max(ylim); end
end
for i = 1:4, 
    subplot(2,2,i); 
    ylim([0 maxylim]); 
%     fill([greybar_start greybar_end greybar_end greybar_start],[min(ylim) min(ylim) max(ylim) max(ylim)], [0.75 0.75 0.75], 'FaceAlpha',0.5,'LineStyle','none');
%     fill([yellowbar_start yellowbar_end yellowbar_end yellowbar_start],[min(ylim) min(ylim) max(ylim) max(ylim)], [1 1 0], 'FaceAlpha',0.5,'LineStyle','none');
    line([endperiod endperiod],ylim,'Color','k','LineStyle',':'); line([startperiod startperiod],ylim,'Color','k','LineStyle',':');
    set(gca, 'Tickdir','Out','xlim',[-300 1300],'XColor','k','YColor','k');
end    
            
save_group_prestim_psth = 0;
if save_group_prestim_psth
    fid = fopen('Z:\LabTools\Matlab\TEMPO_Analysis\ProtocolSpecific\CuedDirectionDiscrim\group_prestimPSTH_shortdurdelay.txt','W');
    fprintf(fid,templabels);
    fprintf(fid,'\r\n');
    for i = 1:size(temp,1)
        fprintf(fid, '%9.5f\t', temp(i,:));
        fprintf(fid, '\r\n');
    end
    fclose(fid)
end

%% sacc psth: this figure contains saccade psths for cued trials
handl(2) = figure;
set(handl(2),'PaperPosition', [.2 .2 8 10.7], 'Position', [250 50 500 573], 'Name', 'Group data: Raw Sacc PSTH');
for i = 1:length(unique_coherence)
    for j = 1:2
        subplot(length(unique_coherence),2, 2*(i-1)+j); hold on
        for k = 1:length(unique_cue_dir_type)
            plot(sacc_x, mean_sacc_psth{i,j,k}, linetypes{k}, 'LineWidth',1.2);
        end
        axis('tight'); upyl(j) = max(ylim); ylim([0 upyl(j)]);
        if j == 1
            ylabel(sprintf('Coh = %3.1f%%',unique_coherence(i)));
        end
        if i == 1
            if j==1
                title ('Pref Choices');
            else
                title ('Null Choices');
            end
        end
            if i == length(unique_coherence)
                xlabel ('Time about Saccade onset (ms)');
            end
    end
    [maxyl, maxylind] = max(upyl);
    for j=1:2
        subplot(length(unique_coherence),2,2*(i-1)+j);
        if j~=maxylind
            ylim([0 maxyl]);
        end
        if plot_ttests
            if ~sum(isnan(h_sacc{i,j}))
                c1 = sacc_x(logical(h_sacc{i,j} & (mean_sacc_psth{i,j,comp(1)} > mean_sacc_psth{i,j,comp(2)})));
                c2 = sacc_x(logical(h_sacc{i,j} & (mean_sacc_psth{i,j,comp(1)} < mean_sacc_psth{i,j,comp(2)})));
                plot(c1, max(ylim).*ones(length(c1),1),sprintf('%s+',comp_color{1}),'MarkerSize',3);
                plot(c2, max(ylim).*ones(length(c2),1),sprintf('%s+',comp_color{2}),'MarkerSize',3);
            end
        end
    end
end
%this is another version of the above, this time averaging across non-0 coherences and saving that
temp = sacc_x';
ilabels = {'0','n0'};
jlabels = {'T1','T2'};
klabels = {'NLc','NUc','PRc'};
templabels = 'Time ';
for j = 1:2
    for i = 1:2
        for k = 1:3
            if i==1 
                temp = [temp mean_sacc_psth{i,j,k}'];
            else %for non zero coherences, average across the non-zero psths and append that to temp
                crossmean_psth = mean([mean_sacc_psth{2,j,k}; mean_sacc_psth{3,j,k}; ...
                      mean_sacc_psth{4,j,k}; mean_sacc_psth{5,j,k}], 1); 
                temp = [temp crossmean_psth'];
            end
            templabels = [templabels, ilabels{i}, jlabels{j}, klabels{k}, ' '];
        end
    end
end
save_group_sacc_psth = 0;
if save_group_sacc_psth
    fid = fopen('Z:\LabTools\Matlab\TEMPO_Analysis\ProtocolSpecific\CuedDirectionDiscrim\group_saccPSTH_shortdurdelay.txt','W');
    fprintf(fid,templabels);
    fprintf(fid,'\r\n');
    for i = 1:size(temp,1)
        fprintf(fid, '%9.5f\t', temp(i,:));
        fprintf(fid, '\r\n');
    end
    fclose(fid)
end


%% Plot the cue / cue-only responses
handl(3) = figure; %this figure contains the combined cueonset response (from motion trials), 
                   %and all three periods from cueonly trials
set(handl(3),'PaperPosition', [.2 .2 8 10.7], 'Position', [250 50 500 573], 'Name', 'Group data: Raw Cue Response PSTH');

%postcue_psth combined across both directions - all cued trials
subplot(5,2,1:2); hold on;
for k = 1:length(unique_cue_dir_type)
    cueh(k) = plot(postcue_x, mean_postcue_combined_psth{k},linetypes{k}, 'LineWidth',1.2);
end
axis('tight'); ylim([0 max(ylim)]); ylim([0 30]);
xlabel('Time about cue onset (ms)'); ylabel('FR (Hz)');
if plot_ttests
    if ~sum(isnan(h_postcue_combined))
        c1 = postcue_x(logical(h_postcue_combined & (mean_postcue_combined_psth{comp(1)} > mean_postcue_combined_psth{comp(2)})));
        c2 = postcue_x(logical(h_postcue_combined & (mean_postcue_combined_psth{comp(1)} < mean_postcue_combined_psth{comp(2)})));
        plot(c1, max(ylim).*ones(length(c1),1),sprintf('%s+',comp_color{1}),'MarkerSize',3);
        plot(c2, max(ylim).*ones(length(c2),1),sprintf('%s+',comp_color{2}),'MarkerSize',3);
    end
end
% temp = [postcue_x' mean_postcue_combined_psth{1}' mean_postcue_combined_psth{2}' mean_postcue_combined_psth{3}'];
% save('Z:\LabTools\Matlab\TEMPO_Analysis\ProtocolSpecific\CuedDirectionDiscrim\group_cueonset_shortdurdelay.txt', '-ascii', 'temp')

%postcue_psth - all cued trials
for j = 1:2
    subplot(5,2,2+j); hold on;
    for k = 1:length(unique_cue_dir_type)
        plot(postcue_x, mean_postcue_psth{j,k},linetypes{k}, 'LineWidth',1.2);
    end
    axis('tight'); upyl(j) = max(ylim); ylim([0 upyl(j)]);
    xlabel('Time about cue onset');
    if j == 1;
        title('Pref Choices: Motion');
        ylabel('FR (Hz)');
    else
        title('Null Choices: Motion');
    end
end
[maxyl, maxylind] = max(upyl);
for j=1:2
    subplot(5,2,2+j);
    if j~=maxylind
        ylim([0 maxyl]);
    end
    if plot_ttests
        if ~sum(isnan(h_postcue{j}))
            c1 = postcue_x(logical(h_postcue{j} & (mean_postcue_psth{j,comp(1)} > mean_postcue_psth{j,comp(2)})));
            c2 = postcue_x(logical(h_postcue{j} & (mean_postcue_psth{j,comp(1)} < mean_postcue_psth{j,comp(2)})));
            plot(c1, max(ylim).*ones(length(c1),1),sprintf('%s+',comp_color{1}),'MarkerSize',3);
            plot(c2, max(ylim).*ones(length(c2),1),sprintf('%s+',comp_color{2}),'MarkerSize',3);
        end
    end
end

%postcue_psth - cueonly trials
for j = 1:2
    subplot(5,2,4+j); hold on;
    for k = 1:2:length(unique_cue_dir_type)
        plot(postcue_co_x, mean_co_postcue_psth{j,k},linetypes{k}, 'LineWidth',1.2);
    end
    axis('tight'); upyl(j) = max(ylim); ylim([0 upyl(j)]);
    xlabel('Time about cue onset');
    if j == 1;
        title('Pref Choices: CueOnly');
        ylabel('FR (Hz)');
    else
        title('Null Choices: CueOnly');
    end
end
[maxyl, maxylind] = max(upyl);
for j=1:2
    subplot(5,2,2+j);
    if j~=maxylind
        ylim([0 maxyl]);
    end
    if plot_ttests
        if ~sum(isnan(h_co_postcue{j}))
            c1 = postcue_co_x(logical(h_co_postcue{j} & (mean_co_postcue_psth{j,comp_co(1)} > mean_co_postcue_psth{j,comp_co(2)})));
            c2 = postcue_co_x(logical(h_co_postcue{j} & (mean_co_postcue_psth{j,comp_co(1)} < mean_co_postcue_psth{j,comp_co(2)})));
            plot(c1, max(ylim).*ones(length(c1),1),sprintf('%s+',comp_color_co{1}),'MarkerSize',3);
            plot(c2, max(ylim).*ones(length(c2),1),sprintf('%s+',comp_color_co{2}),'MarkerSize',3);
        end
    end
end

%prestim_psth - cueonly trials
for j = 1:2
    subplot(5,2,6+j); hold on;
    for k = 1:2:length(unique_cue_dir_type)
        plot(prestim_co_x, mean_co_prestim_psth{j,k},linetypes{k}, 'LineWidth',1.2);
    end
    axis('tight'); upyl(j) = max(ylim); ylim([0 upyl(j)]);
    xlabel('Time about stim onset');
    if j == 1;
        ylabel('FR (Hz)');
    end
end
[maxyl, maxylind] = max(upyl);
for j=1:2
    subplot(5,2,6+j);
    if j~=maxylind
        ylim([0 maxyl]);
    end
    if plot_ttests
        if ~sum(isnan(h_co_prestim{j}))
            c1 = prestim_co_x(logical(h_co_prestim{j} & (mean_co_prestim_psth{j,comp_co(1)} > mean_co_prestim_psth{j,comp_co(2)})));
            c2 = prestim_co_x(logical(h_co_prestim{j} & (mean_co_prestim_psth{j,comp_co(1)} < mean_co_prestim_psth{j,comp_co(2)})));
            plot(c1, max(ylim).*ones(length(c1),1),sprintf('%s+',comp_color_co{1}),'MarkerSize',3);
            plot(c2, max(ylim).*ones(length(c2),1),sprintf('%s+',comp_color_co{2}),'MarkerSize',3);
        end
    end
end

%sacc_psth - cueonly trials
for j = 1:2
    subplot(5,2,8+j); hold on;
    for k = 1:2:length(unique_cue_dir_type)
        plot(sacc_x, mean_co_sacc_psth{j,k},linetypes{k}, 'LineWidth',1.2);
    end
    axis('tight'); upyl(j) = max(ylim); ylim([0 upyl(j)]);
    xlabel('Time about sacc onset');
    if j == 1;
        ylabel('FR (Hz)');
    end
end
[maxyl, maxylind] = max(upyl);
for j=1:2
    subplot(5,2,8+j);
    if j~=maxylind
        ylim([0 maxyl]);
    end
    if plot_ttests
        if ~sum(isnan(h_co_sacc{j}))
            c1 = sacc_x(logical(h_co_sacc{j} & (mean_co_sacc_psth{j,comp_co(1)} > mean_co_sacc_psth{j,comp_co(2)})));
            c2 = sacc_x(logical(h_co_sacc{j} & (mean_co_sacc_psth{j,comp_co(1)} < mean_co_sacc_psth{j,comp_co(2)})));
            plot(c1, max(ylim).*ones(length(c1),1),sprintf('%s+',comp_color_co{1}),'MarkerSize',3);
            plot(c2, max(ylim).*ones(length(c2),1),sprintf('%s+',comp_color_co{2}),'MarkerSize',3);
        end
    end
end

%% also plot the saccadic activity broken down by correct/incorrect
handl(5) = figure; %this figure contains saccade psths for cued trials
set(handl(5),'PaperPosition', [.2 .2 8 10.7], 'Position', [250 50 500 573], 'Name', 'Group data: Raw Sacc PSTH - Correct Sort');
for i = 1:length(unique_coherence)
    for j = 1:2 %choice direction
        for g = 1:2 %motion direction
                
            subplot(length(unique_coherence),4, 4*(i-1) + (g==j)*g + (g~=j)*(j+2)); hold on
            for k = 1:length(unique_cue_dir_type)
                plot(sacc_x, mean_ungrouped_sacc_psth{i,j,k,g}, linetypes{k}, 'LineWidth',1.2);
            end
            axis('tight'); upyl(2*(j-1)+g) = max(ylim); %ylim([0 upyl(j)]);
            if (j==1)&(g==1)
                ylabel(sprintf('Coh = %3.1f%%',unique_coherence(i)));
            end
            if i == 1
                if (j==1)&(g==1), title('Corr Pref Choices'); 
                elseif (j==2)&(g==2), title('Corr Null Choices');
                elseif (j==1)&(g==2), title('Inc Pref Choices');
                else, title('Inc Null Choices');
                end
            end
            if i == length(unique_coherence)
                xlabel ('Time about Saccade onset');
            end
        end
    end
    [maxyl, maxylind] = max(upyl);
    for j=1:4
        subplot(length(unique_coherence),4,4*(i-1)+j);
        ylim([0 maxyl]);
        if 0 %plot_ttests
            if ~sum(isnan(h_sacc{i,j}))
                c1 = sacc_x(logical(h_sacc{i,j} & (mean_sacc_psth{i,j,comp(1)} > mean_sacc_psth{i,j,comp(2)})));
                c2 = sacc_x(logical(h_sacc{i,j} & (mean_sacc_psth{i,j,comp(1)} < mean_sacc_psth{i,j,comp(2)})));
                plot(c1, max(ylim).*ones(length(c1),1),sprintf('%s+',comp_color{1}),'MarkerSize',3);
                plot(c2, max(ylim).*ones(length(c2),1),sprintf('%s+',comp_color{2}),'MarkerSize',3);
            end
        end
    end
end

%% this figure contains motion psths for cued trials, split by motion direction and choice
handl(6) = figure; 
set(handl(6),'PaperPosition', [.2 .2 8 10.7], 'Position', [250 50 500 573], 'Name', 'Group data: Raw Motion PSTH, Motion+Choice split');
for i = 1:length(unique_coherence)
    for j = 1:2 %choice
        for g = 1:2 %prefdir motion=1, nulldir mot=2
            column = 2*(j-1)+g; %T1Ch+T1Mot, T1Ch+T2Mot, T2Ch+T1Mot, T2Ch+T2Mot
            subplot(length(unique_coherence),4,4*(i-1)+column); hold on;
            for k = 1:length(unique_cue_dir_type)
                plot(prestim_x, mean_ungrouped_prestim_psth{i,j,k,g}, linetypes{k}, 'LineWidth',1.2);
            end
            axis('tight'); upyl(column) = max(ylim); ylim([0 upyl(column)]);
            if column == 1
                ylabel(sprintf('Coh = %3.1f%%',unique_coherence(i)));
            end
            if i == 1
                switch column
                    case 1, title('T1 ch + T1 mot');
                    case 2, title('T1 ch + T2 mot');
                    case 3, title('T2 ch + T1 mot');
                    case 4, title('T2 ch + T2 mot');
                end
            end
            if i == length(unique_coherence)
                xlabel ('Time about Motion onset');
            end
        end
    end
    [maxyl, maxylind] = max(upyl);
    for column=1:4
        subplot(length(unique_coherence),4,4*(i-1)+column); hold on
        j = ceil(column/2); g = 2-mod(column,2);
        ylim([0 maxyl]);
        if plot_ttests
            if ~sum(isnan(h_prestim_ungrouped{i,j,g}))
                c1 = prestim_x(logical(h_prestim_ungrouped{i,j,g} & (mean_ungrouped_prestim_psth{i,j,comp(1),g} > mean_ungrouped_prestim_psth{i,j,comp(2),g})));
                c2 = prestim_x(logical(h_prestim_ungrouped{i,j,g} & (mean_ungrouped_prestim_psth{i,j,comp(1),g} < mean_ungrouped_prestim_psth{i,j,comp(2),g})));
                plot(c1, max(ylim).*ones(length(c1),1),sprintf('%s+',comp_color{1}),'MarkerSize',3);
                plot(c2, max(ylim).*ones(length(c2),1),sprintf('%s+',comp_color{2}),'MarkerSize',3);
            end
        end
    end
end
%this is another version of the above, this time averaging across non-0 coherences and saving that
temp = prestim_x'; clear crossmean_psth;
ilabels = {'0','n0'};
jlabels = {'T1','T2'};
klabels = {'NLc','NUc','PRc'};
glabels = {'Pm','Nm'};
templabels = 'Time ';
for j = 1:2
    for g = 1:2
        for i = 1:2
            for k = 1:3
                if i==1
                    temp = [temp mean_ungrouped_prestim_psth{i,j,k,g}'];
                else %for non zero coherences, average across the non-zero psths and append that to temp
                    crossmean_psth{j,k,g} = mean([mean_ungrouped_prestim_psth{2,j,k,g}; mean_ungrouped_prestim_psth{3,j,k,g}; ...
                        mean_ungrouped_prestim_psth{4,j,k,g}; mean_ungrouped_prestim_psth{5,j,k,g}], 1);
                    temp = [temp crossmean_psth{j,k,g}'];
                end
                templabels = [templabels, ilabels{i}, jlabels{j}, klabels{k}, glabels{g}, ' '];
            end
        end
    end
end
figure;
for i = 1:2
    for j = 1:2 %choice
        for g = 1:2 %prefdir motion=1, nulldir mot=2
            column = 2*(j-1)+g; %T1Ch+T1Mot, T1Ch+T2Mot, T2Ch+T1Mot, T2Ch+T2Mot
            subplot(2,4,4*(i-1)+column); hold on;
            fill([greybar_start greybar_end greybar_end greybar_start],[0 0 50 50], [0.75 0.75 0.75], 'FaceAlpha',0.5,'LineStyle','none');
            fill([yellowbar_start yellowbar_end yellowbar_end yellowbar_start],[0 0 50 50], [1 1 0], 'FaceAlpha',0.5,'LineStyle','none');
            line([endperiod endperiod],[0 50],'Color','k','LineStyle',':');
            line([startperiod startperiod],[0 50],'Color','k','LineStyle',':');
            for k = 1:length(unique_cue_dir_type)
                if i==1, plot(prestim_x, mean_ungrouped_prestim_psth{i,j,k,g}, linetypes{k}, 'LineWidth',2);
                else  plot(prestim_x, crossmean_psth{j,k,g}, linetypes{k}, 'LineWidth',2);
                end
            end
            axis('tight'); upyl(column) = max(ylim); ylim([0 upyl(column)]);
        end
    end
    [maxyl, maxylind] = max(upyl);
    for column=1:4
        subplot(2,4,4*(i-1)+column); 
        ylim([0 maxyl]); xlim([min(prestim_x),max(prestim_x)]);
    end
end
save_group_ungrouped_prestim_psth = 0;
if save_group_ungrouped_prestim_psth
    fid = fopen('Z:\LabTools\Matlab\TEMPO_Analysis\ProtocolSpecific\CuedDirectionDiscrim\group_prestimPSTHxMotion_shortdurdelay.txt','W');
    fprintf(fid,templabels);
    fprintf(fid,'\r\n');
    for i = 1:size(temp,1)
        fprintf(fid, '%9.5f\t', temp(i,:));
        fprintf(fid, '\r\n');
    end
    fclose(fid)
end

%% first this figure contains psths for cued trials sorted only by motion and coherence, aligned to motion onset
handl(8) = figure; 
set(handl(8),'PaperPosition', [.2 .2 8 10.7], 'Position', [250 50 500 573], 'Name', 'Group data: Raw Motion PSTH - aligned to Motion and Coher');
for i = 1:length(unique_coherence)
    for g = 1:2
        subplot(length(unique_coherence),2, 2*(i-1)+g); hold on
        for k = 1:length(unique_cue_dir_type)
            plot(prestim_x, mean_motgrouped_prestim_psth{i,k,g}, linetypes{k}, 'LineWidth',1.2);
        end
        axis('tight'); upyl(g) = max(ylim); ylim([0 upyl(g)]);
        if g == 1
            ylabel(sprintf('Coh = %3.1f%%',unique_coherence(i)));
        end
        if i == 1
            if g==1
                title ('Pref Motion','Color','k');
            else
                title ('Null Motion','Color','k');
            end
        end
        if i == length(unique_coherence)
            xlabel ('Time about Motion onset');
        end
    end
    [maxyl, maxylind] = max(upyl);
    for g=1:2
        subplot(length(unique_coherence),2,2*(i-1)+g);
        if g~=maxylind
            ylim([0 maxyl]);
        end
%         if plot_ttests
%             if ~sum(isnan(h_prestim{i,j}))
%                 c1 = prestim_x(logical(h_prestim{i,j} & (mean_prestim_psth{i,j,comp(1)} > mean_prestim_psth{i,j,comp(2)})));
%                 c2 = prestim_x(logical(h_prestim{i,j} & (mean_prestim_psth{i,j,comp(1)} < mean_prestim_psth{i,j,comp(2)})));
%                 plot(c1, max(ylim).*ones(length(c1),1),sprintf('%s+',comp_color{1}),'MarkerSize',3);
%                 plot(c2, max(ylim).*ones(length(c2),1),sprintf('%s+',comp_color{2}),'MarkerSize',3);
%             end
%         end
        set(gca,'tickdir','out','xcolor','k','ycolor','k');
        fill([greybar_start greybar_end greybar_end greybar_start],[min(ylim) min(ylim) max(ylim) max(ylim)], [0.75 0.75 0.75], 'FaceAlpha',0.5,'LineStyle','none');
        fill([yellowbar_start yellowbar_end yellowbar_end yellowbar_start],[min(ylim) min(ylim) max(ylim) max(ylim)], [1 1 0], 'FaceAlpha',0.5,'LineStyle','none');
        line([endperiod endperiod],ylim,'Color','k','LineStyle',':');
        line([startperiod startperiod],ylim,'Color','k','LineStyle',':');
    end
end

%this formats labels and data for easy import into origin.  note that this
%does NOT import the t-test marker locations.
% temp = prestim_x';
% ilabels = {'0','2','4','8','16'};
% glabels = {'T1','T2'};
% klabels = {'NLc','NUc','PRc'};
% templabels = 'Time ';
% for g = 1:2
%     for i = 1:5
%         for k = 1:3
%             temp = [temp mean_motgrouped_prestim_psth{i,k,g}'];
%             templabels = [templabels, ilabels{i}, glabels{g}, klabels{k}, ' '];
%         end
%     end
% end
% fid = fopen('Z:\LabTools\Matlab\TEMPO_Analysis\ProtocolSpecific\CuedDirectionDiscrim\group_motgroupedprestimPSTH_shortdurdelay.txt','W');
% fprintf(fid,templabels);
% fprintf(fid,'\r\n');
% for i = 1:size(temp,1)
%     fprintf(fid, '%9.5f\t', temp(i,:));
%     fprintf(fid, '\r\n');
% end
% fclose(fid)

% this is another version of the above, this time averaging across non-0 coherences and saving that
temp = prestim_x'; clear crossmean_psth;
ilabels = {'0','n0'};
glabels = {'T1','T2'};
klabels = {'NLc','NUc','PRc'};
templabels = 'Time ';
for g = 1:2
    for i = 1:2
        for k = 1:3
            if i==1 
                temp = [temp mean_motgrouped_prestim_psth{i,k,g}'];
            else %for non zero coherences, average across the non-zero psths and append that to temp
                crossmean_psth{g,k} = mean([mean_motgrouped_prestim_psth{2,k,g}; mean_motgrouped_prestim_psth{3,k,g}; ...
                      mean_motgrouped_prestim_psth{4,k,g}; mean_motgrouped_prestim_psth{5,k,g}], 1); 
                temp = [temp crossmean_psth{g,k}'];
            end
            templabels = [templabels, ilabels{i}, glabels{g}, klabels{k}, ' '];
        end
    end
end
figure; hold on;
maxylim = 0;
for g = 1:2
    subplot(2,2,g); hold on; %zero-coh
    fill([greybar_start greybar_end greybar_end greybar_start],[0 0 50 50], [0.75 0.75 0.75], 'FaceAlpha',0.5,'LineStyle','none');
    fill([yellowbar_start yellowbar_end yellowbar_end yellowbar_start],[0 0 50 50], [1 1 0], 'FaceAlpha',0.5,'LineStyle','none');
    for k = 1:3, plot(prestim_x, mean_motgrouped_prestim_psth{1,k,g}, linetypes{k}, 'LineWidth',2); end
    if g==1, ylabel('Firing Rate (Hz)','Color','k'); end
    if max(ylim) > maxylim, maxylim = max(ylim); end
    subplot(2,2,2+g); hold on; %nonzero-coh
    fill([greybar_start greybar_end greybar_end greybar_start],[0 0 50 50], [0.75 0.75 0.75], 'FaceAlpha',0.5,'LineStyle','none');
    fill([yellowbar_start yellowbar_end yellowbar_end yellowbar_start],[0 0 50 50], [1 1 0], 'FaceAlpha',0.5,'LineStyle','none');
    for k = 1:3, plot(prestim_x, crossmean_psth{g,k}, linetypes{k}, 'LineWidth',2); end
    if g==1, ylabel('Firing Rate (Hz)','Color','k'); end
    xlabel('Time about Motion Onset (ms)','Color','k');
    if max(ylim) > maxylim, maxylim = max(ylim); end
end
for i = 1:4, 
    subplot(2,2,i); 
    ylim([0 maxylim]); 
%     fill([greybar_start greybar_end greybar_end greybar_start],[min(ylim) min(ylim) max(ylim) max(ylim)], [0.75 0.75 0.75], 'FaceAlpha',0.5,'LineStyle','none');
%     fill([yellowbar_start yellowbar_end yellowbar_end yellowbar_start],[min(ylim) min(ylim) max(ylim) max(ylim)], [1 1 0], 'FaceAlpha',0.5,'LineStyle','none');
    line([endperiod endperiod],ylim,'Color','k','LineStyle',':'); line([startperiod startperiod],ylim,'Color','k','LineStyle',':');
    set(gca, 'Tickdir','Out','xlim',[-300 1300],'XColor','k','YColor','k');
end    
            
save_group_prestim_psth = 0;
if save_group_prestim_psth
    fid = fopen('Z:\LabTools\Matlab\TEMPO_Analysis\ProtocolSpecific\CuedDirectionDiscrim\group_motgroupedprestimPSTH_shortdurdelay.txt','W');
    fprintf(fid,templabels);
    fprintf(fid,'\r\n');
    for i = 1:size(temp,1)
        fprintf(fid, '%9.5f\t', temp(i,:));
        fprintf(fid, '\r\n');
    end
    fclose(fid)
end

%% this figure contains motion psths for cued trials, sort by RT and ISI
handl(7) = figure; 
set(handl(7),'PaperPosition', [.2 .2 8 10.7], 'Position', [250 50 500 573], 'Name', 'Group data: Motion PSTH split by RT (|=ShortRT, :=LongRT)');
for i = 1:length(unique_coherence)
    for j = 1:2
        subplot(length(unique_coherence),2, 2*(i-1)+j); hold on
        for k = 1:length(unique_cue_dir_type)
            plot(prestim_x, mean_lort_prestim_psth{i,j,k}, linetypesT{1,k}, 'LineWidth',1.2);
            plot(prestim_x, mean_hirt_prestim_psth{i,j,k}, linetypesT{2,k}, 'LineWidth',1.2);
        end
        axis('tight'); upyl(j) = max(ylim); ylim([0 upyl(j)]);
        if j == 1
            ylabel(sprintf('Coh = %3.1f%%',unique_coherence(i)));
        end
        if i == 1
            if j==1
                title ('Pref Choices');
            else
                title ('Null Choices');
            end
        end
        if i == length(unique_coherence)
            xlabel ('Time about Motion onset');
        end
    end
    [maxyl, maxylind] = max(upyl);
    for j=1:2
        subplot(length(unique_coherence),2,2*(i-1)+j);
        if j~=maxylind
            ylim([0 maxyl]);
        end
        if plot_ttests && 0
            if ~sum(isnan(h_prestim{i,j}))
                c1 = prestim_x(logical(h_prestim{i,j} & (mean_prestim_psth{i,j,comp(1)} > mean_prestim_psth{i,j,comp(2)})));
                c2 = prestim_x(logical(h_prestim{i,j} & (mean_prestim_psth{i,j,comp(1)} < mean_prestim_psth{i,j,comp(2)})));
                plot(c1, max(ylim).*ones(length(c1),1),sprintf('%s+',comp_color{1}),'MarkerSize',3);
                plot(c2, max(ylim).*ones(length(c2),1),sprintf('%s+',comp_color{2}),'MarkerSize',3);
            end
        end
    end
end

handl(8) = figure; %this figure contains motion psths for cued trials
set(handl(8),'PaperPosition', [.2 .2 8 10.7], 'Position', [250 50 500 573], 'Name', 'Group data: Motion PSTH split by ISI (|=ShortISI, :=LongISI)');
for i = 1:length(unique_coherence)
    for j = 1:2
        subplot(length(unique_coherence),2, 2*(i-1)+j); hold on
        for k = 1:length(unique_cue_dir_type)
            plot(prestim_x, mean_loisi_prestim_psth{i,j,k}, linetypesT{1,k}, 'LineWidth',1.2);
            plot(prestim_x, mean_hiisi_prestim_psth{i,j,k}, linetypesT{2,k}, 'LineWidth',1.2);
        end
        axis('tight'); upyl(j) = max(ylim); ylim([0 upyl(j)]);
        if j == 1
            ylabel(sprintf('Coh = %3.1f%%',unique_coherence(i)));
        end
        if i == 1
            if j==1
                title ('Pref Choices');
            else
                title ('Null Choices');
            end
        end
        if i == length(unique_coherence)
            xlabel ('Time about Motion onset');
        end
    end
    [maxyl, maxylind] = max(upyl);
    for j=1:2
        subplot(length(unique_coherence),2,2*(i-1)+j);
        if j~=maxylind
            ylim([0 maxyl]);
        end
        if plot_ttests && 0
            if ~sum(isnan(h_prestim{i,j}))
                c1 = prestim_x(logical(h_prestim{i,j} & (mean_prestim_psth{i,j,comp(1)} > mean_prestim_psth{i,j,comp(2)})));
                c2 = prestim_x(logical(h_prestim{i,j} & (mean_prestim_psth{i,j,comp(1)} < mean_prestim_psth{i,j,comp(2)})));
                plot(c1, max(ylim).*ones(length(c1),1),sprintf('%s+',comp_color{1}),'MarkerSize',3);
                plot(c2, max(ylim).*ones(length(c2),1),sprintf('%s+',comp_color{2}),'MarkerSize',3);
            end
        end
    end
end

%% now compute  the average differential activity in the last part of the
% motion onset on the preferred choice trials...

begin_ramp = 700; end_ramp = 1100;
inds = find( (prestim_x >= begin_ramp-1) & (prestim_x <= end_ramp) );
figure;
for i = 1:length(unique_coherence)
    subplot(length(unique_coherence),1,i); hold on;
    for k = 1:3
        ramp{i,k} = diff(mean_prestim_psth{i,1,k}(inds));
        plot(begin_ramp:end_ramp, ramp{i,k}, linetypes{k}, 'LineWidth',1.2);
    end
    p_anova_ramp(i,:) = anova2([diff(cum_prestim_psth{i,1,1}(gooddata,:),1,2); diff(cum_prestim_psth{i,1,2}(gooddata,:),1,2); ...
    diff(cum_prestim_psth{i,1,3}(gooddata,:),1,2)], n, 'off');
%     if (i==3)
%         ylabel(sprintf('Differential firing rate (Hz/s)\nCoh=4%%'))
%     else
    ylabel(sprintf('Coh=%d%%',unique_coherence(i)));
%     end
end



%% now compute the latency - use whatever smoothing i have already
% implemented and 
offset = 0; %HACK to make sure cutoff isn't too early
inds = find( prestim_x > offset );
preinds = find(prestim_x <= offset );
goodcells = find(gooddata);
for i = 1:length(unique_coherence)
    for k = 1:3
        pref_psth = cum_prestim_psth{i,1,k}(gooddata,:);
        null_psth = cum_prestim_psth{i,2,k}(gooddata,:);
        diff_psth = pref_psth - null_psth; 
        diff_baseline = diff_psth(:,preinds);
        diff_thresh = mean(diff_baseline,2) + 3.*std(diff_baseline,0,2);
        for m = 1:sum(gooddata)
            temp = find( (diff_psth(gooddata(m),inds) > diff_thresh(m)) ,1);
            if isempty(temp)
                lat{i,k}(m) = max(inds);
            else
                lat{i,k}(m) = temp + offset;
            end
        end        
        diff_meanpsth = mean_prestim_psth{i,1,k} - mean_prestim_psth{i,2,k};
        diff_meanthresh = mean(diff_meanpsth(preinds)) + 3.*std(diff_meanpsth(preinds));
        mean_lat(i,k) = find( diff_meanpsth(inds) > diff_meanthresh, 1) + offset;
    end
end
offset = 0;
figure
for i = 1:length(unique_coherence)
    subplot(5,1,i); hold on
    for k = 1:3
        diff_meanpsth = mean_prestim_psth{i,1,k} - mean_prestim_psth{i,2,k};
        t = 0; keep_looping = 1;
        while keep_looping
            preinds = find(prestim_x < t);
            thresh = mean(diff_meanpsth(preinds)) + 2.*std(diff_meanpsth(preinds));
            if ( (diff_meanpsth(prestim+t) > thresh) | (t == max(prestim_x)) );
                keep_looping = 0;
                mean_lat(i,k) = t + offset;
            else
                t = t+1;
            end
        end        
        plot(prestim_x,diff_meanpsth,linetypes{k}, 'LineWidth',1.2);
    end
    axis tight;
    plot(repmat(mean_lat(i,1),2,1),ylim,linetypes{1});
    plot(repmat(mean_lat(i,2),2,1),ylim,linetypes{2});
    plot(repmat(mean_lat(i,3),2,1),ylim,linetypes{3});
end



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %now plot the normalized data
% 
% handl(4) = figure;
% set(handl(4),'PaperPosition', [.2 .2 8 10.7], 'Position', [250 50 500 573], 'Name', 'Group data: Norm Motion PSTH');
% for i = 1:length(unique_coherence)
%     for j = 1:2
%         subplot(length(unique_coherence),2, 2*(i-1)+j); hold on
%         for k = 1:length(unique_cue_dir_type)
%             plot(prestim_x, normmean_prestim_psth{i,j,k}, linetypes{k}, 'LineWidth',1.2);
%         end
%         axis('tight'); upyl(j) = max(ylim); ylim([0 upyl(j)]);
%         if j == 1
%             ylabel(sprintf('Coh = %3.1f%%',unique_coherence(i)));
%         end
%         if i == 1
%             if j==1
%                 title ('Pref Choices');
%             else
%                 title ('Null Choices');
%             end
%         end
%         if i == length(unique_coherence)
%             xlabel ('Time about Motion onset');
%         end
%     end
%     [maxyl, maxylind] = max(upyl);
%     for j=1:2
%         subplot(length(unique_coherence),2,2*(i-1)+j);
%         if j~=maxylind
%             ylim([0 maxyl]);
%         end
%         if plot_ttests
%             if ~sum(isnan(h_norm_prestim{i,j}))
%                 plot(prestim_x(logical(h_norm_prestim{i,j})), max(ylim).*ones(sum(h_norm_prestim{i,j})),'m+','MarkerFaceColor',[1 1 0],'MarkerSize',3);
%             end
%         end
%     end
% end
% 
% handl(5) = figure;
% set(handl(5),'PaperPosition', [.2 .2 8 10.7], 'Position', [250 50 500 573], 'Name', 'Group data: Norm Sacc PSTH');
% for i = 1:length(unique_coherence)
%     for j = 1:2
%         subplot(length(unique_coherence),2, 2*(i-1)+j); hold on
%         for k = 1:length(unique_cue_dir_type)
%             plot(sacc_x, normmean_sacc_psth{i,j,k}, linetypes{k}, 'LineWidth',1.2);
%         end
%         axis('tight'); upyl(j) = max(ylim); ylim([0 upyl(j)]);
%         if j == 1
%             ylabel(sprintf('Coh = %3.1f%%',unique_coherence(i)));
%         end
%         if i == 1
%             if j==1
%                 title ('Pref Choices');
%             else
%                 title ('Null Choices');
%             end
%         end
%         if i == length(unique_coherence)
%             xlabel ('Time about Saccade onset');
%         end
%     end
%     [maxyl, maxylind] = max(upyl);
%     for j=1:2
%         subplot(length(unique_coherence),2,2*(i-1)+j);
%         if j~=maxylind
%             ylim([0 maxyl]);
%         end
%         if plot_ttests
%             if ~sum(isnan(h_norm_sacc{i,j}))
%                 plot(sacc_x(logical(h_norm_sacc{i,j})), max(ylim).*ones(sum(h_norm_sacc{i,j})),'m+','MarkerFaceColor',[1 1 0],'MarkerSize',3);
%             end
%         end
%     end
% end
% 
% handl(6) = figure;
% set(handl(6),'PaperPosition', [.2 .2 8 10.7], 'Position', [250 50 500 573], 'Name', 'Group data: Norm Cue Response PSTH');
% 
% %postcue_psth combined across both directions - all cued trials
% subplot(5,2,1:2); hold on;
% for k = 1:length(unique_cue_dir_type)
%     plot(postcue_x, normmean_postcue_combined_psth{k},linetypes{k}, 'LineWidth',1.2);
% end
% axis('tight'); ylim([0 max(ylim)]);
% xlabel('Time about cue onset'); ylabel('FR (Hz)');
% if plot_ttests
%     if ~sum(isnan(h_norm_postcue_combined))
%         plot(postcue_x(logical(h_norm_postcue_combined)), max(ylim).*ones(sum(h_norm_postcue_combined)),'m+','MarkerFaceColor',[1 1 0],'MarkerSize',3);
%     end
% end
% 
% %postcue_psth - all cued trials
% for j = 1:2
%     subplot(5,2,2+j); hold on;
%     for k = 1:length(unique_cue_dir_type)
%         plot(postcue_x, normmean_postcue_psth{j,k},linetypes{k}, 'LineWidth',1.2);
%     end
%     axis('tight'); upyl(j) = max(ylim); ylim([0 upyl(j)]);
%     xlabel('Time about cue onset');
%     if j == 1;
%         title('Pref Choices: Motion');
%         ylabel('FR (Hz)');
%     else
%         title('Null Choices: Motion');
%     end
% end
% [maxyl, maxylind] = max(upyl);
% for j=1:2
%     subplot(5,2,2+j);
%     if j~=maxylind
%         ylim([0 maxyl]);
%     end
%     if plot_ttests
%         if ~sum(isnan(h_norm_postcue{j}))
%             plot(postcue_x(logical(h_norm_postcue{j})), max(ylim).*ones(sum(h_norm_postcue{j})),'m+','MarkerFaceColor',[1 1 0],'MarkerSize',3);
%         end
%     end
% end
% 
% %postcue_psth - cueonly trials
% for j = 1:2
%     subplot(5,2,4+j); hold on;
%     for k = 1:2:length(unique_cue_dir_type)
%         plot(postcue_x, normmean_co_postcue_psth{j,k},linetypes{k}, 'LineWidth',1.2);
%     end
%     axis('tight'); upyl(j) = max(ylim); ylim([0 upyl(j)]);
%     xlabel('Time about cue onset');
%     if j == 1;
%         title('Pref Choices: CueOnly');
%         ylabel('FR (Hz)');
%     else
%         title('Null Choices: CueOnly');
%     end
% end
% [maxyl, maxylind] = max(upyl);
% for j=1:2
%     subplot(5,2,2+j);
%     if j~=maxylind
%         ylim([0 maxyl]);
%     end
%     if plot_ttests
%         if ~sum(isnan(h_norm_co_postcue{j}))
%             plot(postcue_x(logical(h_norm_co_postcue{j})), max(ylim).*ones(sum(h_norm_co_postcue{j})),'m+','MarkerFaceColor',[1 1 0],'MarkerSize',3);
%         end
%     end
% end
% 
% %prestim_psth - cueonly trials
% for j = 1:2
%     subplot(5,2,6+j); hold on;
%     for k = 1:2:length(unique_cue_dir_type)
%         plot(prestim_x, normmean_co_prestim_psth{j,k},linetypes{k}, 'LineWidth',1.2);
%     end
%     axis('tight'); upyl(j) = max(ylim); ylim([0 upyl(j)]);
%     xlabel('Time about stim onset');
%     if j == 1;
%         ylabel('FR (Hz)');
%     end
% end
% [maxyl, maxylind] = max(upyl);
% for j=1:2
%     subplot(5,2,6+j);
%     if j~=maxylind
%         ylim([0 maxyl]);
%     end
%     if plot_ttests
%         if ~sum(isnan(h_norm_co_prestim{j}))
%             plot(prestim_x(logical(h_norm_co_prestim{j})), max(ylim).*ones(sum(h_norm_co_prestim{j})),'m+','MarkerFaceColor',[1 1 0],'MarkerSize',3);
%         end
%     end
% end
% 
% %sacc_psth - cueonly trials
% for j = 1:2
%     subplot(5,2,8+j); hold on;
%     for k = 1:2:length(unique_cue_dir_type)
%         plot(sacc_x, normmean_co_sacc_psth{j,k},linetypes{k}, 'LineWidth',1.2);
%     end
%     axis('tight'); upyl(j) = max(ylim); ylim([0 upyl(j)]);
%     xlabel('Time about sacc onset');
%     if j == 1;
%         ylabel('FR (Hz)');
%     end
% end
% [maxyl, maxylind] = max(upyl);
% for j=1:2
%     subplot(5,2,8+j);
%     if j~=maxylind
%         ylim([0 maxyl]);
%     end
%     if plot_ttests
%         if ~sum(isnan(h_norm_co_sacc{j}))
%             plot(sacc_x(logical(h_norm_co_sacc{j})), max(ylim).*ones(sum(h_norm_co_sacc{j})),'m+','MarkerFaceColor',[1 1 0],'MarkerSize',3);
%         end
%     end
% end
% 
