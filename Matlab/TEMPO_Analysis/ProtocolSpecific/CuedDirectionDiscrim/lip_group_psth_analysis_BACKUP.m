%this loads all each output file and combines it into 

% % for Baskin's shortdur data use the following
% analysis_dir = 'Z:\Data\Tempo\Baskin\Analysis\LIP_PSTH\ShortDur';
% d = ls(analysis_dir);
% d = d(3:end,:); %removes '.' and '..' from top of list
% exclude_cells = []; %2 as example

% %for Baskin's shortdur+delay data use the following
% analysis_dir = 'Z:\Data\Tempo\Baskin\Analysis\LIP_PSTH\ShortDurDelay';
% d = ls(analysis_dir);
% d = d(3:end,:); %removes '.' and '..' from top of list
% % d = d([2:end],:); %kluge to remove cell 1 which has problems
% exclude_cells = [2 15 18 28]; 

% %for Baskin's full dur(1000ms) data use the following:
% analysis_dir = 'Z:\Data\Tempo\Baskin\Analysis\LIP_PSTH\';
% d = ls(analysis_dir);
% d = d(6:end,:); %removes '.' and '..' and 3 subdirectory names from top of list
% d = d([1:60 62:length(d)],:) %KLUGE to remove cell 61 which has problems 
% exclude_cells = [1 4 9 15 17 31 40 43 47 51 52 53 57 60];

cum_file = []; cum_coher = []; cum_postcue_combined_psth = repmat({[]},[1,3]); cum_normval = [];
cum_postcue_psth = repmat({[]},[2,3]); cum_prestim_psth = repmat({[]},[5,2,3]); cum_sacc_psth = repmat({[]},[5,2,3]); 
cum_co_postcue_psth = repmat({[]},[2,3]); cum_co_prestim_psth = repmat({[]},[2,3]); cum_co_sacc_psth = repmat({[]},[2,3]); 
cum_ungrouped_prestim_psth = repmat({[]},[5,2,3,2]); cum_ungrouped_sacc_psth = repmat({[]},[5,2,3,2]); 
cum_zero_pct_xing = [];
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
% SAVEFILE = 'Z:\LabTools\Matlab\TEMPO_Analysis\ProtocolSpecific\CuedDirectionDiscrim\cum_lip_psth_summary.mat';
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
    cum_file{length(cum_file)+1} = file;
    cum_coher = [cum_coher; coher];
    cum_normval = [cum_normval; normval];
    cum_zero_pct_xing = [cum_zero_pct_xing; offset];
    for k = 1:3 %length(unique_cue_dir_type)
        temp = conv(gaussfilt, long_sm_postcue_combined_psth{k});
        cum_postcue_combined_psth{k}(end+1,:) = temp(lb+buff+1:end-lb-buff);
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
                end
            end
        end
    end
end





% DATAFILE = 'Z:\LabTools\Matlab\TEMPO_Analysis\ProtocolSpecific\CuedDirectionDiscrim\cum_lip_psth_TEST.mat'

% load(DATAFILE);
n_total = size(cum_postcue_psth{1,1},1);
gooddata = logical(ones(n_total,1)); %these indicies indicate the subset of datasets that are well-isolated and show decent behavior
gooddata(exclude_cells) = 0;
gooddata_co = gooddata;
for i = 1:length(gooddata_co) %exclude cells where the monkey didn't err at least once on both directions on the cueonly trials
    if ( isnan(cum_co_postcue_psth{1,1}(i,1)) || isnan(cum_co_postcue_psth{2,1}(i,1)) || ...
            isnan(cum_co_postcue_psth{2,3}(i,1)) || isnan(cum_co_postcue_psth{2,3}(i,1)) ) 
        gooddata_co(i) = 0;
    end
end
n = sum(gooddata);
n_co = sum(gooddata_co);


unique_cue_dir_type = [-1 0 1]; %-1 = NullCue, 0 = Neutral Cue, 1 = PrefCue
unique_cue_dir_type_names = {'NullCue','Neutral','PrefCue'};
unique_coherence = cum_coher(1,:);
unique_direction = [1 2]; %1=PrefDir motion, 2=NullDir motion
%these timings are determined in PSTH_CuedDirec
precue = 200; %time to show before the cue starts
postcue = 400; %time to show after cue starts
prestim = 300; %time to display before the visual stimulus
poststim = 1100; %time to display after the visual stimulus
presacc = 400; %time to show before the saccade starts
postsacc = 200; %time to show after saccade starts
precue_co = 200;
postcue_co = 150+1000;
prestim_co = -800; %200 before FP off
poststim_co = 1400;
postcue_x = [-precue:(size(cum_postcue_psth{1,1,1},2)-precue-1)];
prestim_x = [-prestim:(size(cum_prestim_psth{1,1,1},2)-prestim-1)];
postcue_co_x = [-precue_co:(length(cum_co_postcue_psth{1,1})-precue_co-1)];
prestim_co_x = [-200:(length(cum_co_prestim_psth{1,1})-200-1)]; %relative to FP offset

sacc_x = [-presacc:(size(cum_sacc_psth{1,1},2)-presacc-1)];
linetypes = {'b-','r-','g-'};

clear temp cued_sacc_onset co_sacc_onset combined_sacc_onset;
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

%first compute means and normalize means for all psth.  Normalization is
%to the target onset of the corresponding trial type
%also for those data sets where the monkey didn't make both possible
%choices for the cue only trials.  
for k = 1:length(unique_cue_dir_type)
    
    select = ~isnan(cum_postcue_combined_psth{k}(:,1));
    mean_postcue_combined_psth{k} = mean(cum_postcue_combined_psth{k}(gooddata & select,:),1);
    normval_postcue_combined_psth{k} = repmat(normval,1,size(cum_postcue_combined_psth{k},2));
    norm_postcue_combined_psth{k} = cum_postcue_combined_psth{k}./normval_postcue_combined_psth{k};
    normmean_postcue_combined_psth{k} = mean(norm_postcue_combined_psth{k}(gooddata & select,:),1);

    for j = 1:2
        select = ~isnan(cum_postcue_psth{j,k}(:,1));
        mean_postcue_psth{j,k} = mean(cum_postcue_psth{j,k}(gooddata & select,:),1);
        normval_postcue_psth{j,k} = repmat(normval,1,size(cum_postcue_psth{j,k},2));
        norm_postcue_psth{j,k} = cum_postcue_psth{j,k}./normval_postcue_psth{j,k};
        normmean_postcue_psth{j,k} = mean(norm_postcue_psth{j,k}(gooddata & select,:),1);

        if k ~= 2
            select = ~isnan(cum_co_postcue_psth{j,k}(:,1));
            mean_co_postcue_psth{j,k} = mean(cum_co_postcue_psth{j,k}(gooddata_co & select,:),1);
            normval_co_postcue_psth{j,k} = repmat(normval,1,size(cum_co_postcue_psth{j,k},2));
            norm_co_postcue_psth{j,k} = cum_co_postcue_psth{j,k}./normval_co_postcue_psth{j,k};
            normmean_co_postcue_psth{j,k} = mean(norm_co_postcue_psth{j,k}(gooddata_co & select,:),1);

            select = ~isnan(cum_co_prestim_psth{j,k}(:,1));
            mean_co_prestim_psth{j,k} = mean(cum_co_prestim_psth{j,k}(gooddata_co & select,:),1);
            normval_co_prestim_psth{j,k} = repmat(normval,1,size(cum_co_prestim_psth{j,k},2));
            norm_co_prestim_psth{j,k} = cum_co_prestim_psth{j,k}./normval_co_prestim_psth{j,k};
            normmean_co_prestim_psth{j,k} = mean(norm_co_prestim_psth{j,k}(gooddata_co & select,:),1);

            select = ~isnan(cum_co_sacc_psth{j,k}(:,1));
            mean_co_sacc_psth{j,k} = mean(cum_co_sacc_psth{j,k}(gooddata_co & select,:),1);
            normval_co_sacc_psth{j,k} = repmat(normval,1,size(cum_co_sacc_psth{j,k},2));
            norm_co_sacc_psth{j,k} = cum_co_sacc_psth{j,k}./normval_co_sacc_psth{j,k};
            normmean_co_sacc_psth{j,k} = mean(norm_co_sacc_psth{j,k}(gooddata_co & select,:),1);
        end

        for i = 1:length(unique_coherence)
            select = ~isnan(cum_prestim_psth{i,j,k}(:,1));
            mean_prestim_psth{i,j,k} = mean(cum_prestim_psth{i,j,k}(gooddata & select,:),1);
            normval_prestim_psth{i,j,k} = repmat(normval,1,size(cum_prestim_psth{i,j,k},2));
            norm_prestim_psth{i,j,k} = cum_prestim_psth{i,j,k}./normval_prestim_psth{i,j,k};
            normmean_prestim_psth{i,j,k} = mean(norm_prestim_psth{i,j,k}(gooddata & select,:),1);

            select = ~isnan(cum_sacc_psth{i,j,k}(:,1));
            mean_sacc_psth{i,j,k} = mean(cum_sacc_psth{i,j,k}(gooddata & select,:),1);
            normval_sacc_psth{i,j,k} = repmat(normval,1,size(cum_sacc_psth{i,j,k},2));
            norm_sacc_psth{i,j,k} = cum_sacc_psth{i,j,k}./normval_sacc_psth{i,j,k};
            normmean_sacc_psth{i,j,k} = mean(norm_sacc_psth{i,j,k}(gooddata & select,:),1);
            
            for g = 1:2
                temp = gooddata;
                temp(excluded_cells_ungrouped_data{i,j,k,g}) = 0; %exclude those cells that lacked the corresponding condition
                
                mean_ungrouped_prestim_psth{i,j,k,g} = mean(cum_ungrouped_prestim_psth{i,j,k,g}(temp,:),1);
                normval_ungrouped_prestim_psth{i,j,k,g} = repmat(normval,1,size(cum_ungrouped_prestim_psth{i,j,k,g},2));
                norm_ungrouped_prestim_psth{i,j,k,g} = cum_ungrouped_prestim_psth{i,j,k,g}./normval_ungrouped_prestim_psth{i,j,k,g};
                normmean_ungrouped_prestim_psth{i,j,k,g} = mean(norm_ungrouped_prestim_psth{i,j,k,g}(temp,:),1);

                mean_ungrouped_sacc_psth{i,j,k,g} = mean(cum_ungrouped_sacc_psth{i,j,k,g}(temp,:),1);
                normval_ungrouped_sacc_psth{i,j,k,g} = repmat(normval,1,size(cum_ungrouped_sacc_psth{i,j,k,g},2));
                norm_ungrouped_sacc_psth{i,j,k,g} = cum_ungrouped_sacc_psth{i,j,k,g}./normval_ungrouped_sacc_psth{i,j,k,g};
                normmean_ungrouped_sacc_psth{i,j,k,g} = mean(norm_ungrouped_sacc_psth{i,j,k,g}(temp,:),1);
            end
        end
    end
end

%setup data(y) and covariate(x) matrices for the ramp and baserate portions
%of the cell-by-cell and population data.  This prepares 4 regressions.
begin_ramp = 700; end_ramp = 1000;
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
                slope = regress(cum_ungrouped_prestim_psth{i,1,k,g}(temp(m),rampinds)',[ones(length(rampinds),1) (1:length(rampinds))']);
                slope_neu = regress(cum_ungrouped_prestim_psth{i,1,2,g}(temp(m),rampinds)',[ones(length(rampinds),1) (1:length(rampinds))']);
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


%I'm not sure but i believe the following code is to orgnaize the data
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

%test for changes in offset for motion-aligned PSTHs
%first subtract neutral cue response for each cell. 
%second, average across all conditions (5 coherences * 2 choices)
%NB: matlab's ttest ignores those elements in the array with value NaN (so don't need to worry about cells that don't have a trials in a particular condition)
counter = 0; clear meanpsthdiff grouppsthdiff psthdiff;
t_on = -100; t_off = 0;
for i = 1:length(unique_coherence)
    for j = 1:2 %choice
        counter = counter + 1;
        for k = [1 3] %only operate on T1- and T2-cued trials
            psthdiff{i,j,(k+1)/2} = cum_prestim_psth{i,j,k}(gooddata,prestim+t_on+1:prestim+t_off) - cum_prestim_psth{i,j,2}(gooddata,prestim+t_on+1:prestim+t_off);
            meanpsthdiff{i,j}(:,(k+1)/2) = mean(psthdiff{i,j,(k+1)/2},2); %averages across time in psth into delta firing rates (with funny units)
            grouppsthdiff{(k+1)/2}(:,(i-1)*2+j) = meanpsthdiff{i,j}(:,(k+1)/2); 
            %reorganize the average firing rates so that each matrix has a
            %list of the firing rate differences for the 10 combinations of coherence and choice.  This is for subsequent testing.
        end
        [h_offset(i,j) p_offset(i,j)] = ttest(meanpsthdiff{i,j}(:,1), meanpsthdiff{i,j}(:,2));
    end
end
[h_offset p_offset]
[h_groupoffset p_groupoffset] = ttest(mean(grouppsthdiff{1},2), mean(grouppsthdiff{2},2))
[nanmean(nanmean(grouppsthdiff{1},2)) nanmean(nanmean(grouppsthdiff{2},2))]

%test for changes in slope for motion-aligned PSTHs
t_on = 700; t_off = 1000;
for i = 1:length(unique_coherence)
    counter = 0;
    for m = 1:length(cum_file)
        if gooddata(m)
            counter = counter+1;
            for k = 1:3
                temp = cum_prestim_psth{i,1,k}(m,prestim+t_on+1:prestim+t_off)';
                if isnan(temp(1))
                    cue_slope{i,k}(counter) = NaN;
                else
                    b_temp = regress(temp, [ones(size(temp)) [1:t_off-t_on]']);
                    cue_slope{i,k}(counter) = b_temp(2);
                end
            end
        end
    end
    diff_slope{i}(:,1) = cue_slope{i,1}-cue_slope{i,2};
    diff_slope{i}(:,2) = cue_slope{i,3}-cue_slope{i,2};
    [h_slope(i) p_slope(i)] = ttest(diff_slope{i}(:,1), diff_slope{i}(:,2));
    groupdiff_slope{1}(i,:) = diff_slope{i}(:,1)'; %within each matrix, each row has different cells' slope for one coherence
    groupdiff_slope{2}(i,:) = diff_slope{i}(:,2)';
end
[h_slope' p_slope']
meandiff_slope = [nanmean(groupdiff_slope{1},1)' nanmean(groupdiff_slope{2},1)'];
[h_groupslope p_groupslope] = ttest(meandiff_slope(:,1), meandiff_slope(:,2))
mean(meandiff_slope,1)

%compute CMI (Cue Modulation Index) for Valid and Invalid trials for each cell and plot scatter
%CMI = (cue - neutral) / (cue + neutral), taking the mean response during the specified window
t_on = 0; t_off = 250;
clear CMI meanCMI;
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
[h,p] = ttest(meanCMI(1,:),meanCMI(2,:))
[mean(meanCMI,2) std(meanCMI,0,2)]
figure;
plot(meanCMI(1,:),meanCMI(2,:),'b*')
xlabel('T2-cue CMI'); ylabel('T1-cue CMI');
title(sprintf('CMI: %dms to %dms post-motion onset', t_on, t_off));
%xlim([-.2,.2]); ylim([-.2,.2]);
[h,p] = ttest(meanCMI(2,:),0)
[h,p] = ttest(meanCMI(1,:),0)
% keyboard;

                
% figure %plot estimates of latency
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

%plot neutral cue psths on same axes at different coherences
h = figure ;
set(h,'PaperPosition', [.2 .2 8 10.7], 'Position', [250 50 500 573], 'Name', 'PSTH by Coherence');
coh_lines = {'r','m','b','c','g'};
coh_dashedlines = {'r--','m--','b--','c--','g--'};
for k = 1:3
    legstr = 'legh(k) = legend(coh_h';
    subplot(3,1,k); hold on;
    for i = 1:length(unique_coherence)
%         coh_h(i) = plot(prestim_x, mean_ungrouped_prestim_psth{i,1,k,1},coh_lines{i});
        coh_h(i) = plot(prestim_x, mean_prestim_psth{i,1,k},coh_lines{i},'LineWidth',1.2);
        plot(prestim_x, mean_prestim_psth{i,2,k},coh_dashedlines{i},'LineWidth',1.2);
        legstr = strcat(legstr, sprintf(',''%d%%''',unique_coherence(i)));
    end
    axis tight;
    legstr = strcat(legstr, ',''Location'',''NorthWest'');');
    eval(legstr); set(legh(k),'box','off');
    title(sprintf('PrefChoice Trials: %s',unique_cue_dir_type_names{k}));
    ylabel('Firing Rate (Hz)');
    if k==3 
        xlabel('Time about Motion Onset'); 
    end
end
%now test whether slopes truly are different prior to T1-choices in neutral trials
lo_bnd = prestim+1+700; hi_bnd = prestim+1100;
counter = 0;
for m = 1:length(cum_file)
    if gooddata(m)
        counter = counter+1;
        for i = 1:length(unique_coherence)
            temp = cum_prestim_psth{i,1,2}(m,lo_bnd:hi_bnd); 
            b_temp = regress(temp', [ones(size(temp')) [lo_bnd-prestim:hi_bnd-prestim]']);
            neuslope(counter,i) = b_temp(2);
            neuoffset(counter,i) = b_temp(1);
        end
    end
end

% keyboard 
% %%%%%%%%%%%%%
% %now perform paired t-tests comparing prefdir cues and neutral cues
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

%now do permutation test, shuffling across assignment of two cuedirs in comp
nboot = 1000; alpha = 0.05;
CIlobnd = floor(nboot*(alpha/2)); CIupbnd = ceil(nboot*(1-alpha/2)); %demarcate central 95%
for b = 1:nboot
    goodcells = find(gooddata);
    for i = 1:length(goodcells)
        temp = ceil(2*rand);
        boot_cum_psth{1}(i,:) = cum_prestim_psth{comp(temp)}(goodcells(i),:);
        boot_cum_psth{2}(i,:) = cum_prestim_psth{comp(3-temp)}(goodcells(i),:);
    end
    boot_psth{1}(b,:) = nanmean(boot_cum_psth{1},1);
    boot_psth{2}(b,:) = nanmean(boot_cum_psth{2},1);
end
for v = 1:length(prestim_x)
    temp = sort(boot_psth{1}(:,v));
    prestim_psth_CI{1}(:,v) = [temp(CIlobnd); temp(CIupbnd)];
    temp = sort(boot_psth{2}(:,v));
    prestim_psth_CI{2}(:,v) = [temp(CIlobnd); temp(CIupbnd)];
end


%anovas to compare activity across timecourses... 
p_anova_postcue_combined = anova2([cum_postcue_combined_psth{1}(gooddata,:); cum_postcue_combined_psth{2}(gooddata,:); cum_postcue_combined_psth{3}(gooddata,:)], n,'off');
for j = 1:2
    p_anova_postcue(j,:) = anova2([cum_postcue_psth{j,1}(gooddata,:); cum_postcue_psth{j,2}(gooddata,:); cum_postcue_psth{j,3}(gooddata,:)], n,'off');
%     p_anova_co_postcue(j,:) = anova2([cum_co_postcue_psth{j,1}(gooddata_co,:); cum_co_postcue_psth{j,3}(gooddata_co,:)], n_co,'off');
%     p_anova_co_prestim(j,:) = anova2([cum_co_prestim_psth{j,1}(gooddata_co,:); cum_co_prestim_psth{j,3}(gooddata_co,:)], n_co,'off');
%     p_anova_co_sacc(j,:) = anova2([cum_co_sacc_psth{j,1}(gooddata_co,:); cum_co_sacc_psth{j,3}(gooddata_co,:)], n_co,'off');
    for i = 1:length(coher)
        p_anova_prestim(i,j,:) = anova2([cum_prestim_psth{i,j,1}(gooddata,:); cum_prestim_psth{i,j,2}(gooddata,:); cum_prestim_psth{i,j,3}(gooddata,:)], n,'off');
        p_anova_sacc(i,j,:) = anova2([cum_sacc_psth{i,j,1}(gooddata,:); cum_sacc_psth{i,j,2}(gooddata,:); cum_sacc_psth{i,j,3}(gooddata,:)], n,'off');
    end
end
%repeated measure anovas (type help rm_anova2 for explanation of variables)
y = [reshape(cum_postcue_combined_psth{1}(gooddata,:)',prod(size(cum_postcue_combined_psth{1}(gooddata,:))),1); reshape(cum_postcue_combined_psth{2}(gooddata,:)',prod(size(cum_postcue_combined_psth{2}(gooddata,:))),1); reshape(cum_postcue_combined_psth{3}(gooddata,:)',prod(size(cum_postcue_combined_psth{3}(gooddata,:))),1)];
s = reshape(repmat([1:n], size(cum_postcue_combined_psth{1},2), 3),size(y));
f1 = repmat([1:size(cum_postcue_combined_psth{1},2)]', length(y)/size(cum_postcue_combined_psth{1},2),1);
f2 = reshape(repmat([1:3], prod(size(cum_postcue_combined_psth{1}(gooddata,:))),1),size(y));
f1name = 'Time';  f2name = 'CueDirType';
% p_rm_anova_postcue_combined_full = rm_anova2(y,s,f1,f2,{f1name, f2name});
% p_rm_anova_postcue_combined = [p_rm_anova_postcue_combined{3,6} p_rm_anova_postcue_combined_full{4,6}]; %[p_cuedirtype p_time_X_cuedirtype]
for j = 1:2
    y = [reshape(cum_postcue_psth{j,1}(gooddata,:)',prod(size(cum_postcue_psth{j,1}(gooddata,:))),1); reshape(cum_postcue_psth{j,2}(gooddata,:)',prod(size(cum_postcue_psth{j,2}(gooddata,:))),1); reshape(cum_postcue_psth{j,3}(gooddata,:)',prod(size(cum_postcue_psth{j,3}(gooddata,:))),1)];
    s = reshape(repmat([1:n], size(cum_postcue_psth{j,1},2), 3),size(y));
    f1 = repmat([1:size(cum_postcue_psth{j,1},2)]', length(y)/size(cum_postcue_psth{j,1},2), 1);
    f2 = reshape(repmat([1:3], prod(size(cum_postcue_psth{j,1}(gooddata,:))),1),size(y));
%     p_rm_anova_postcue_full{j} = rm_anova2(y,s,f1,f2,{f1name, f2name});
%     p_rm_anova_postcue(j,:) = [p_rm_anova_postcue_full{j}{3,6} p_rm_anova_postcue_full{j}{4,6}]
    for i = 1:length(coher)
        y = [reshape(cum_prestim_psth{i,j,1}(gooddata,:)',prod(size(cum_prestim_psth{i,j,1}(gooddata,:))),1); reshape(cum_prestim_psth{i,j,2}(gooddata,:)',prod(size(cum_prestim_psth{i,j,2}(gooddata,:))),1); reshape(cum_prestim_psth{i,j,3}(gooddata,:)',prod(size(cum_prestim_psth{i,j,3}(gooddata,:))),1)];
        s = reshape(repmat([1:n], size(cum_prestim_psth{i,j,1},2), 3),size(y));
        f1 = repmat([1:size(cum_prestim_psth{i,j,1},2)]', length(y)/size(cum_prestim_psth{i,j,1},2), 1);
        f2 = reshape(repmat([1:3], prod(size(cum_prestim_psth{i,j,1}(gooddata,:))),1),size(y));
%         p_rm_anova_prestim_full{i,j} = rm_anova2(y,s,f1,f2,{f1name, f2name});
%         p_rm_anova_prestim(i,j,:) = [p_rm_anova_prestim_full{i,j}{3,6} p_rm_anova_prestim_full{i,j}{4,6}]
        
        y = [reshape(cum_sacc_psth{i,j,1}(gooddata,:)',prod(size(cum_sacc_psth{i,j,1}(gooddata,:))),1); reshape(cum_sacc_psth{i,j,2}(gooddata,:)',prod(size(cum_sacc_psth{i,j,2}(gooddata,:))),1); reshape(cum_sacc_psth{i,j,3}(gooddata,:)',prod(size(cum_sacc_psth{i,j,3}(gooddata,:))),1)];
        s = reshape(repmat([1:n], size(cum_sacc_psth{i,j,1},2), 3),size(y));
        f1 = repmat([1:size(cum_sacc_psth{i,j,1},2)]', length(y)/size(cum_sacc_psth{i,j,1},2), 1);
        f2 = reshape(repmat([1:3], prod(size(cum_sacc_psth{i,j,1}(gooddata,:))),1),size(y));
%         p_rm_anova_sacc_full{i,j} = rm_anova2(y,s,f1,f2,{f1name, f2name});
%         p_rm_anova_sacc(i,j,:) = [p_rm_anova_sacc_full{i,j}{3,6} p_rm_anova_sacc_full{i,j}{4,6}]
    end
	y = [reshape(cum_co_postcue_psth{j,1}(gooddata_co,:)',prod(size(cum_co_postcue_psth{j,1}(gooddata_co,:))),1); reshape(cum_co_postcue_psth{j,3}(gooddata_co,:)',prod(size(cum_co_postcue_psth{j,3}(gooddata_co,:))),1)];
    s = reshape(repmat([1:n_co], size(cum_co_postcue_psth{1},2), 2),size(y));
    f1 = repmat([1:size(cum_co_postcue_psth{1},2)]', length(y)/size(cum_postcue_psth{1},2),1);
    f2 = reshape(repmat([1 2], prod(size(cum_co_postcue_psth{1}(gooddata_co,:))),1),size(y));
%     p_rm_anova_co_postcue_full{j} = rm_anova2(y,s,f1,f2,{f1name, f2name});
%     p_rm_anova_co_postcue(j,:) = [p_rm_anova_co_postcue_full{j}{3,6} p_rm_anova_co_postcue_full{j}{4,6}]
    
	y = [reshape(cum_co_prestim_psth{j,1}(gooddata_co,:)',prod(size(cum_co_prestim_psth{j,1}(gooddata_co,:))),1); reshape(cum_co_prestim_psth{j,3}(gooddata_co,:)',prod(size(cum_co_prestim_psth{j,3}(gooddata_co,:))),1)];
    s = reshape(repmat([1:n_co], size(cum_co_prestim_psth{1},2), 2),size(y));
    f1 = repmat([1:size(cum_co_prestim_psth{1},2)]', length(y)/size(cum_prestim_psth{1},2),1);
    f2 = reshape(repmat([1 2], prod(size(cum_co_prestim_psth{1}(gooddata_co,:))),1),size(y));
%     p_rm_anova_co_prestim_full{j} = rm_anova2(y,s,f1,f2,{f1name, f2name});
%     p_rm_anova_co_prestim(j,:) = [p_rm_anova_co_prestim_full{j}{3,6} p_rm_anova_co_prestim_full{j}{4,6}]
    
	y = [reshape(cum_co_sacc_psth{j,1}(gooddata_co,:)',prod(size(cum_co_sacc_psth{j,1}(gooddata_co,:))),1); reshape(cum_co_sacc_psth{j,3}(gooddata_co,:)',prod(size(cum_co_sacc_psth{j,3}(gooddata_co,:))),1)];
    s = reshape(repmat([1:n_co], size(cum_co_sacc_psth{1},2), 2),size(y));
    f1 = repmat([1:size(cum_co_sacc_psth{1},2)]', length(y)/size(cum_sacc_psth{1},2),1);
    f2 = reshape(repmat([1 2], prod(size(cum_co_sacc_psth{1}(gooddata_co,:))),1),size(y));
%     p_rm_anova_co_sacc_full{j} = rm_anova2(y,s,f1,f2,{f1name, f2name});
%     p_rm_anova_co_sacc(j,:) = [p_rm_anova_co_sacc_full{j}{3,6} p_rm_anova_co_sacc_full{j}{4,6}]
end    

% keyboard

plot_ttests = 1;

%%%%%%%%%%%%%
%now plot these averaged raw data
handl(1) = figure; %this figure contains motion psths for cued trials
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
        if plot_ttests
            if ~sum(isnan(h_prestim{i,j}))
                c1 = prestim_x(logical(h_prestim{i,j} & (mean_prestim_psth{i,j,comp(1)} > mean_prestim_psth{i,j,comp(2)})));
                c2 = prestim_x(logical(h_prestim{i,j} & (mean_prestim_psth{i,j,comp(1)} < mean_prestim_psth{i,j,comp(2)})));
                plot(c1, max(ylim).*ones(length(c1),1),sprintf('%s+',comp_color{1}),'MarkerSize',3);
                plot(c2, max(ylim).*ones(length(c2),1),sprintf('%s+',comp_color{2}),'MarkerSize',3);
            end
        end
    end
end

handl(2) = figure; %this figure contains saccade psths for cued trials
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

%also plot the saccadic activity broken down by correct/incorrect
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

%%%%%%%%%%%%%%%
% now compute  the average differential activity in the last part of the
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



%%%%%%%%%%%%
% now compute the latency - use whatever smoothing i have already
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
