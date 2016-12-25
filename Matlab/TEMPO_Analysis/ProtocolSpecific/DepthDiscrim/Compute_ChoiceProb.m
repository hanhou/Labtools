%-----------------------------------------------------------------------------------------------------------------------
%-- Compute_ChoiceProb.m -- Uses ROC analysis to compute a choice probability for each different stimulus level
%--	GCD, 5/30/00
%-----------------------------------------------------------------------------------------------------------------------
function [grandCP, grandPval] = Compute_ChoiceProb(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

TEMPO_Defs;		%needed for defines like IN_T1_WIN_CD
Path_Defs;
ProtocolDefs;	%needed for all protocol specific functions - contains keywords - BJP 1/4/01

%disp('computing choice probabilities...');

Pref_HDisp = data.one_time_params(PREFERRED_HDISP);

%get the column of values of horiz. disparities in the dots_params matrix
h_disp = data.dots_params(DOTS_HDISP,:,PATCH1);
unique_hdisp = munique(h_disp');

%get the binocular correlations
binoc_corr = data.dots_params(DOTS_BIN_CORR, :, PATCH1);
unique_bin_corr = munique(binoc_corr');

%get signed binocular correlations
sign = (h_disp == Pref_HDisp)*2 - 1;	%=1 if preferred disparity, -1 if null disparity
signed_bin_corr = binoc_corr .* sign;
unique_signed_bin_corr = munique(signed_bin_corr');
%[h_disp' sign' binoc_corr' signed_bin_corr']

%get the random seed for each trial of the Patch1 dots
%check to see if there is a fixed seed and store this for later if there is.
if (size(data.dots_params,1) >= DOTS_BIN_CORR_SEED)  %for backwards compatibility with old files that lack this
    seeds = data.dots_params(DOTS_BIN_CORR_SEED, :, PATCH1);
    select_fixed_seeds = logical(seeds == data.one_time_params(FIXED_SEED));
else 
    select_fixed_seeds = [];
end
if (sum(select_fixed_seeds) >= 1)
    fixed_seed = data.one_time_params(FIXED_SEED);
else
    fixed_seed = NaN;
end

%now, get the firing rates for all the trials 
spike_rates = data.spike_rates(SpikeChan, :);
%start_offset = -200; % start of calculation relative to stim onset, ms
%window_size = 200;  % window size, ms
%spike_rates = ComputeSpikeRates(data, length(h_disp), StartCode, StartCode, start_offset+30, start_offset+window_size+30);
%spike_rates = spike_rates(1,:);
    
%get indices of any NULL conditions (for measuring spontaneous activity
null_trials = logical( (binoc_corr == data.one_time_params(NULL_VALUE)) );

%now, select trials that fall between BegTrial and EndTrial
trials = 1:length(binoc_corr);		% a vector of trial indices
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );

%now, determine the choice that was made for each trial, PREFERRED or NULL
%by definition, a preferred choice will be made to Target1 and a null choice to Target 2
%thus, look for the events IN_T1_WIN_CD and IN_T2_WIN_CD.  GCD, 5/30/2000
num_trials = length(binoc_corr);
PREFERRED = 1;
NULL = 2;
for i=1:num_trials
    temp = data.event_data(1,:,i);
    events = temp(temp>0);  % all non-zero entries
    if (sum(events == IN_T1_WIN_CD) > 0)
        choice(i) = PREFERRED;
    elseif (sum(events == IN_T2_WIN_CD) > 0)
        choice(i) = NULL;
    else
        disp('Neither T1 or T2 chosen.  This should not happen!.  File must be bogus.');
    end        
end

CORRECT_FOR_VERGENCE = 0;
if (CORRECT_FOR_VERGENCE)    
    %now, Z-score the spike rates for each bin_corr and disparity condition
    %These Z-scored responses will be used to remove the effects of vergence angle
    Z_Spikes = spike_rates;
    for i=1:length(unique_bin_corr)
        for j=1:length(unique_hdisp)
            select = (binoc_corr == unique_bin_corr(i)) & (h_disp == unique_hdisp(j));
            z_dist = spike_rates(select);
            z_dist = (z_dist - mean(z_dist))/std(z_dist);
            Z_Spikes(select) = z_dist;
        end
    end
    %now, get the vergence data for each trial and regress this against Z_Spikes
    [h_verg, v_verg, h_conj, v_conj, calib_h_verg, calib_v_verg, calib_h_conj, calib_v_conj] = ...
        DDiscrim_GetEyeData(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    
    X_fit = [ones(length(calib_h_verg'),1) calib_h_verg'];
    [b, bint, r, rint, stats] = regress(Z_Spikes', X_fit);
    stats;
    Verg_Zspike_P= stats(3);
    
    figure;
    plot(calib_h_verg', Z_Spikes', 'ro');
    hold on;
    plot(calib_h_verg', r, 'go');
    spike_rates = r';
end

CORRECT_FOR_SLOW_SPIKE_RATE_CHANGE = 0;
if (CORRECT_FOR_SLOW_SPIKE_RATE_CHANGE)    
    %now, Z-score the spike rates for each bin_corr and disparity condition
    %These Z-scored responses will be used to remove the effects of slow spike rate change
    Z_Spikes = spike_rates;
    for i=1:length(unique_bin_corr)
        for j=1:length(unique_hdisp)
            select = (binoc_corr == unique_bin_corr(i)) & (h_disp == unique_hdisp(j));
            z_dist = spike_rates(select);
            z_dist = (z_dist - mean(z_dist))/std(z_dist);
            Z_Spikes(select) = z_dist;
        end
    end
    
    %Do a regression of Zspikes against trial number.
    trial_temp = [trials; ones(1,length(trials))];
    [b, bint, r, rint, stats] = regress(Z_Spikes', trial_temp');
    
    figure;
    set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [50 120 500 573], 'Name', 'Zspikes');
    subplot(2, 1, 1);
    hold on;
    Handl(1) = plot(trials(choice == PREFERRED)', Z_Spikes(choice == PREFERRED)', 'ko', 'MarkerFaceColor', 'k');
    hold on;
    Handl(2) = plot(trials(choice == NULL)', Z_Spikes(choice == NULL)', 'ko');
    xlabel('Trials');
    ylabel('Z-Scores');
    legend(Handl, 'Preferred', 'Null', 2);
    titl = sprintf('File: %s', FILE);
    title(titl);
    
    subplot(2, 1, 2);
    hold on;
    Handl(1) = plot(trials(choice == PREFERRED)', r(choice == PREFERRED)', 'ko', 'MarkerFaceColor', 'k');
    hold on;
    Handl(2) = plot(trials(choice == NULL)', r(choice == NULL)', 'ko');
    xlabel('Trials');
    ylabel('Z-Scores');
    legend(Handl, 'Preferred', 'Null', 2);
    
    spike_rates = r';
end

%now, plot the spike distributions, sorted by choice, for each correlation level
figure;
set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [50 120 500 573], 'Name', 'Choice Probabilities');
num_corrs = length(unique_bin_corr);
choice_prob = [];

for i=1:num_corrs
    for j=1:length(unique_hdisp)
        subplot(num_corrs, length(unique_hdisp), (i-1)*length(unique_hdisp) + j);
        pref_choices = ( (choice == PREFERRED) & (binoc_corr == unique_bin_corr(i)) & (h_disp == unique_hdisp(j)) );
        pref_dist{i,j} = spike_rates(pref_choices & select_trials);
        null_choices = ( (choice == NULL) & (binoc_corr == unique_bin_corr(i)) & (h_disp == unique_hdisp(j)) );
        null_dist{i,j} = spike_rates(null_choices & select_trials);
        
        %plot the distributions.  This uses a function (in CommonTools) that I wrote.  GCD
        PlotTwoHists(pref_dist{i,j}, null_dist{i,j});
        
        if (i==1)
            ttl = sprintf('disparity = %5.3f', unique_hdisp(j) );
            title(ttl);
        end
        if (j==1)
            lbl = sprintf('%5.1f %%', unique_bin_corr(i) );
            ylabel(lbl);
        end
        
        if ( (length(pref_dist{i,j}) > 0) & (length(null_dist{i,j}) > 0) )
            [choice_prob(i, j), choice_prob_Pval(i, j)] = ROC_signif_test(pref_dist{i,j}, null_dist{i,j});
            cp = sprintf('%5.2f', choice_prob(i,j));
            xl = XLim; yl = YLim;
            text(xl(2), yl(2)/2, cp);
        end
    end
end    

str = sprintf('%s: PrDisp=%5.3f  FixedSeed=%d', FILE, Pref_HDisp, fixed_seed );
xlabel(str);

%pref_dist{i,j} and null_dist{i,j} are cell arrays that hold the preferred and null choice
%distributions for each correlation level and each disparity.
%NOW, we want to Z-score the distributions (preferred and null choices together) and combine across
%correlations and/or disparities.  GCD, 8/10/00
for i=1:num_corrs
    for j=1:length(unique_hdisp)
        %for each condition, combine the preferred and null choices into one dist., then find mean and std
        all_choices = [];
        all_choices = [pref_dist{i,j}  null_dist{i,j}];
        mean_val(i,j) = mean(all_choices);
        std_val(i,j) = std(all_choices);
        %now use the mean_val and std_val to Z-score the original distributions and store separately
        Z_pref_dist{i,j} = (pref_dist{i,j} - mean_val(i,j))/std_val(i,j);
        Z_null_dist{i,j} = (null_dist{i,j} - mean_val(i,j))/std_val(i,j);
    end
end    

%Now, combine the data across disparity at each correlation value and plot distributions again
figure;
set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [550 120 400 573], 'Name', 'Choice Probabilities combined across Disparity');
for i=1:num_corrs
    subplot(num_corrs, 1, i);
    
    %now, combine z-scored data across disparity at each correlation level
    Zpref{i} = []; Znull{i} = [];
    for j=1:length(unique_hdisp)
        %only include this correlation value in grand if the monkey had less than a 3:1 ratio
        %of choices to the two targets (avoid conditions where there are few errors or bad biases)
        if (min(length(Z_pref_dist{i,j}),length(Z_null_dist{i,j})) / max(length(Z_pref_dist{i,j}),length(Z_null_dist{i,j})) >= (1/3) )
            Zpref{i} = [Zpref{i} Z_pref_dist{i,j}];
            Znull{i} = [Znull{i} Z_null_dist{i,j}];
        end
    end
    
    %plot the distributions.  This uses a function (in CommonTools) that I wrote.  GCD
    PlotTwoHists(Zpref{i}, Znull{i});
    
    lbl = sprintf('%5.1f %%', unique_bin_corr(i) );
    ylabel(lbl);
    
    if ( (length(Zpref{i}) > 0) & (length(Znull{i}) > 0) )
        ch_prob(i) = rocN(Zpref{i}, Znull{i}, 100);
        cp = sprintf('%5.2f', ch_prob(i));
        xl = XLim; yl = YLim;
        text(xl(2), yl(2)/2, cp);    
    end
end

%Now, combine the data across correlation values for each disparity
    
%now, combine z-scored data across disparity at each correlation level
Zpref_pdisp = []; Znull_pdisp = []; Zpref_ndisp = []; Znull_ndisp = [];
for i=1:length(unique_hdisp)
    for j=1:num_corrs
        if (unique_bin_corr(j) ~= 0)
            %only include this correlation value in grand if the monkey had less than a 3:1 ratio
            %of choices to the two targets (avoid conditions where there are few errors or bad biases)
            if (min(length(Z_pref_dist{j,i}),length(Z_null_dist{j,i})) / max(length(Z_pref_dist{j,i}),length(Z_null_dist{j,i})) >= (1/3) )
                if (unique_hdisp(i) == Pref_HDisp)
                    Zpref_pdisp = [Zpref_pdisp Z_pref_dist{j,i}];
                    Znull_pdisp = [Znull_pdisp Z_null_dist{j,i}];
                else
                    Zpref_ndisp = [Zpref_ndisp Z_pref_dist{j,i}];
                    Znull_ndisp = [Znull_ndisp Z_null_dist{j,i}];                   
                end
            end
        end
    end
end    
    
if ( (length(Zpref_pdisp) > 0) & (length(Znull_pdisp) > 0) )
    ch_prob_pdisp = rocN(Zpref_pdisp, Znull_pdisp, 100);
end
if ( (length(Zpref_ndisp) > 0) & (length(Znull_ndisp) > 0) )
    ch_prob_ndisp = rocN(Zpref_ndisp, Znull_ndisp, 100);
end


%get a significance value for the overall CP at zero correlation
zeroCP = NaN; zeroPval = NaN;
if ( (length(Zpref{1}) > 0) & (length(Znull{1}) > 0) )
    [zeroCP, zeroPval] = ROC_signif_test(Zpref{1}, Znull{1});
end

%Now, combine data across correlation to get a grand choice probability, and plot distributions again
figure;
set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [600 20 400 300], 'Name', 'Grand Choice Probability');
Zpref_grand = []; Znull_grand = [];
%combine data across correlations into grand distributions
for i=1:num_corrs 
    Zpref_grand = [Zpref_grand Zpref{i}];   
    Znull_grand = [Znull_grand Znull{i}];   
end
PlotTwoHists(Zpref_grand, Znull_grand);

%do permutation test to get P value for grand CP
[grandCP, grandPval] = ROC_signif_test(Zpref_grand, Znull_grand);
titl = sprintf('grand CP = %5.3f, P = %6.4f', grandCP, grandPval);
title(titl);

% %time course of choice probability
% figure;
% set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [250 50 500 573], 'Name', 'Time Course of Choice Probability');
% subplot(2, 1, 1);
% 
% start_offset = -100; % start of calculation relative to stim onset, ms
% window_size = 100;  % window size, ms
% window_step = 20; % step of sliding window, ms
% start_time = start_offset: window_step: 1500;
% choice_prob_tc = [];
% for j = 1:length(start_time) %calculate spike rates for different window
%     Zpref_tc = []; Znull_tc = [];
%     spike_rates = ComputeSpikeRates(data, length(h_disp), StartCode, StartCode, start_time(j)+30, start_time(j)+window_size+30);
%     spike_rates = spike_rates(1,:);
% 
%     for i=1:length(unique_bin_corr)%loop through each binocular correlation levels, and calculate CPs for each
%         for k=1:length(unique_hdisp)%loop through each disparity level.
%             pref_choices = ( (choice == PREFERRED) & (binoc_corr == unique_bin_corr(i)) & (h_disp == unique_hdisp(k)) );
%             pref_dist_tc{j,i,k} = spike_rates(pref_choices & select_trials);
%             null_choices = ( (choice == NULL) & (binoc_corr == unique_bin_corr(i)) & (h_disp == unique_hdisp(k)) );
%             null_dist_tc{j,i,k} = spike_rates(null_choices & select_trials);
%                
%             if ( (length(pref_dist_tc{j,i,k}) > 0) & (length(null_dist_tc{j,i,k}) > 0) )
%                 choice_prob_tc(j,i,k) = rocN(pref_dist_tc{j,i,k}, null_dist_tc{j,i,k}, 100);
%             end
%             
%             %Z-score using means and variances calculated from the whole 1.5sec visual stimulation period
%             %Z_pref_dist_tc{j,i,k} = (pref_dist_tc{j,i,k} - mean_val(i,k))/std_val(i,k);
%             %Z_null_dist_tc{j,i,k} = (null_dist_tc{j,i,k} - mean_val(i,k))/std_val(i,k);
%             %for each condition, combine the preferred and null choices into one dist., then find mean and std
%             all_choices = [];
%             all_choices = [pref_dist_tc{j,i,k}  null_dist_tc{j,i,k}];
%             mean_val(i,k) = mean(all_choices);
%             std_val(i,k) = std(all_choices);
%             Z_pref_dist_tc{j,i,k} = (pref_dist_tc{j,i,k} - mean_val(i,k))/std_val(i,k);
%             Z_null_dist_tc{j,i,k} = (null_dist_tc{j,i,k} - mean_val(i,k))/std_val(i,k);
%             %only include this correlation value if the monkey had less than a 3:1 ratio
%             %of choices to the two targets (avoid conditions where there are few errors or bad biases)            
%             if (min(length(Z_pref_dist_tc{j,i,k}),length(Z_null_dist_tc{j,i,k})) / max(length(Z_pref_dist_tc{j,i,k}),length(Z_null_dist_tc{j,i,k})) >= (1/3) )
%                 Zpref_tc = [Zpref_tc Z_pref_dist_tc{j,i,k}];
%                 Znull_tc = [Znull_tc Z_null_dist_tc{j,i,k}];
%             end 
%         end
%     end
%     choice_prob_tc_norm(j) = rocN(Zpref_tc, Znull_tc, 100);
% end
% bin_center = start_time + (window_size/2);
% 
% %plot choice probabilities against time for each correlation and disparity level
% for i=1:length(unique_bin_corr)
%     for k=1:length(unique_hdisp)
%         if ( (length(pref_dist_tc{1,i,k}) > 0) & (length(null_dist_tc{1,i,k}) > 0) )
%             hold on;
%             plot(bin_center, choice_prob_tc(:,i,k), 'k-');
%             hold off;               
%         end
%     end
% end
% xlabel('Time after Stimulus Onset (ms)');
% ylabel('Choice Probability');
% YLim([0.0 1.0]);
% XLim([-100 1500]);
% 
% %plot grand choice probabilities against time
% subplot(2, 1, 2);
% hold on;
% plot(bin_center, choice_prob_tc_norm, 'k-');
% hold off;
% xlabel('Time after Stimulus Onset (ms)');
% ylabel('Choice Probability');
% YLim([0.0 1.0]);
% XLim([-100 1500]);
% 
% %time course of firing to preferred and null choice
% figure;
% set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [250 50 500 573], 'Name', 'Time Course of Firing');
% subplot(2, 1, 1);
% 
% start_offset = -100; % start of calculation relative to stim onset, ms
% window_size = 20;  % window size, ms
% window_step = 20; % step of sliding window, ms
% start_time = start_offset: window_step: 1500;
% for j = 1:length(start_time) %calculate spike rates for different window
%     spike_rates = ComputeSpikeRates(data, length(h_disp), StartCode, StartCode, start_time(j)+30, start_time(j)+window_size+30);
%     spike_rates = spike_rates(1,:);
% 
%     for i=1:length(unique_bin_corr)%loop through each binocular correlation levels, and calculate CPs for each
%         for k=1:length(unique_hdisp)%loop through each disparity level.
%             pref_choices = ( (choice == PREFERRED) & (binoc_corr == unique_bin_corr(i)) & (h_disp == unique_hdisp(k)) );
%             null_choices = ( (choice == NULL) & (binoc_corr == unique_bin_corr(i)) & (h_disp == unique_hdisp(k)) );
%             
%             %only include this correlation value if the monkey had less than a 3:1 ratio
%             %of choices to the two targets (avoid conditions where there are few errors or bad biases)      
%             if (min(sum(pref_choices),sum(null_choices)) / max(sum(pref_choices),sum(null_choices)) >= (1/3) )
%                 mean_pref(j,i,k) = mean(spike_rates(pref_choices & select_trials));
%                 mean_null(j,i,k) = mean(spike_rates(null_choices & select_trials));
%             else
%                 mean_pref(j,i,k) = NaN;
%                 mean_null(j,i,k) = NaN;
%             end 
%         end 
%     end
% end
% 
% %normalize firing rates
% max_firing = max([max(mean_pref);max(mean_null)]);
% for j = 1:length(start_time)
%     mean_pref_norm(j,:,:) = mean_pref(j,:,:)./max_firing;
%     mean_null_norm(j,:,:) = mean_null(j,:,:)./max_firing;
% end
%        
% bin_center = start_time + (window_size/2);
% 
% %plot normalize firing rate against time for each correlation and disparity level
% for i=1:length(unique_bin_corr)
%     for k=1:length(unique_hdisp)
%         hold on;
%         plot(bin_center, mean_pref_norm(:,i,k), 'k-');
%         plot(bin_center, mean_null_norm(:,i,k), 'k--');
%         hold off;               
%     end
% end
% xlabel('Time after Stimulus Onset (ms)');
% ylabel('Normalized Response');
% YLim([0.0 1.0]);
% XLim([-100 1500]);
% 
% %average across correlation and disparity values
% corr_count = 0;
% for i=1:length(unique_bin_corr)%loop through each binocular correlation levels, and calculate CPs for each
%     for k=1:length(unique_hdisp)%loop through each disparity level.
%         if (isnan(mean_pref_norm(1,i,k)) == 0)
%             corr_count = corr_count + 1;
%             norm_pref(:,corr_count) = mean_pref_norm(:,i,k);
%             norm_null(:,corr_count) = mean_null_norm(:,i,k);
%         end
%     end 
% end
% mean_norm_pref = mean(norm_pref,2);
% mean_norm_null = mean(norm_null,2);
% 
% %plot average normalized firing rate against time
% subplot(2, 1, 2);
% hold on;
% plot(bin_center, mean_norm_pref, 'k-');
% plot(bin_center, mean_norm_null, 'k-');
% hold off;
% xlabel('Time after Stimulus Onset (ms)');
% ylabel('Normalized Response');
% YLim([0.0 1.0]);
% XLim([-100 1500]);

save ('Z:\LabTools\Matlab\TEMPO_Analysis\ProtocolSpecific\DepthDiscrim\CP_data.mat')
disp ('workspace saved')






%-----------------------------------------------------------------------------------------------------------------------------------------------------------
%now print out some summary parameters to the screen and to a cumulative file
pref_indx = find(unique_hdisp == Pref_HDisp);	%index to preferred disparity (which has no var conditions, if any)
null_indx = find(unique_hdisp ~= Pref_HDisp);	%index to null disparity 
if isnan(fixed_seed)	% this run didn't have NOVAR conditions
    str = sprintf('%s %6.2f %6s %6s %6.4f %7.5f %6.4f %7.5f %6.4f %6.4f', FILE, unique_bin_corr(1), '--', '--', zeroCP, zeroPval, grandCP, grandPval, ch_prob_pdisp, ch_prob_ndisp);      
else
    str = sprintf('%s %6.2f %6.4f %6.4f %6.4f %7.5f %6.4f %7.5f %6.4f %6.4f', FILE, unique_bin_corr(1), choice_prob(1,pref_indx), choice_prob(1,null_indx), zeroCP, zeroPval, grandCP, grandPval, ch_prob_pdisp, ch_prob_ndisp);      
end
disp(str);

% printflag = 0;
% outfile = [BASE_PATH 'ProtocolSpecific\DepthDiscrim\CPSummary.dat'];
% if (exist(outfile, 'file') == 0)    %file does not yet exist
%     printflag = 1;
% end
% fsummid = fopen(outfile, 'a');
% if (printflag)
%     fprintf(fsummid, 'FILE\t lo_corr\t CPnovar\t CPvar\t CPzero\t Pzero\t CPgrnd\t Pgrnd\t CPpref\t CPnull\t');
%     fprintf(fsummid, '\r\n');
% end
% fprintf(fsummid, str);
% fprintf(fsummid, '\r\n');
% fclose(fsummid);
% 
% %----------------------------------------------------------------------------------------------------------------------------------------------------------------
% %now, print out spike rates for all trials at the lowest correlation, sorted by choice
% output = 0;
% if (output)
%     i = size(PATH,2) - 1;
%     while PATH(i) ~='\'	%Analysis directory is one branch below Raw Data Dir
%         i = i - 1;
%     end   
%     PATHOUT = [PATH(1:i) 'Analysis\NeuroPsychoCurves\'];
%     i = size(FILE,2) - 1;
%     while FILE(i) ~='.'
%         i = i - 1;
%     end
%     FILEOUT = [FILE(1:i) 'choice_data'];
%     
%     fileid = [PATHOUT FILEOUT];
%     fwriteid = eval(['fopen(fileid, ''w'')']);
%     
%     pref_choices = ( (choice == PREFERRED) & (binoc_corr == unique_bin_corr(1)) );
%     pref_dist = spike_rates(pref_choices & select_trials);
%     null_choices = ( (choice == NULL) & (binoc_corr == unique_bin_corr(1)) );
%     null_dist = spike_rates(null_choices & select_trials);
%     len = [length(pref_dist) length(null_dist)];
%     max_vals = max(len);
%     min_vals = min(len);
%     
%     fprintf(fwriteid,'Pref\tNull\n');  
%     for i=1:min_vals
%         fprintf(fwriteid, '%6.3f\t%6.3f\n', pref_dist(i), null_dist(i));
%     end
%     for i=min_vals+1:max_vals
%         if (length(pref_dist) == max_vals)
%             fprintf(fwriteid, '%6.3f\t\n', pref_dist(i));
%         else
%             fprintf(fwriteid, '\t%6.3f\n', null_dist(i));
%         end
%     end
%     
%     fclose(fwriteid);
% end
% 
% %-----------------------------------------------------------------------------------------------------------------------------------------------------------
% %now print out data for time course of choice probabilities TU 1/23/01
% output2 = 0;
% if (output2)
%     outfile2 = [BASE_PATH 'ProtocolSpecific\DepthDiscrim\CPtime_course.dat'];
%     printflag = 0;
%     if (exist(outfile2, 'file') == 0)    %file does not yet exist
%         printflag = 1;
%     end
%     fid = fopen(outfile2, 'a');
%     if (printflag)
%         fprintf(fid, 'FILE\t Nthr100\t Nthr200\t Nthr300\t Nthr400\t Nthr500\t Nthr600\t Nthr700\t Nthr800\t Nthr900\t Nthr1000\t Nthr1100\t Nthr1200\t Nthr1300\t Nthr1400\t Nthr1500\t ');
%         fprintf(fid, '\r\n');
%         printflag = 0;
%     end
% 
%     buff = sprintf('%s\t ', FILE);
%     for i=1:length(choice_prob_tc_norm)
%         buff = sprintf('%s\t %6.4f\t', buff, choice_prob_tc_norm(i));
%     end
%     fprintf(fid, '%s', buff);
%     fprintf(fid, '\r\n');
%     fclose(fid);
% end
% 
% %-----------------------------------------------------------------------------------------------------------------------------------------------------------
% %now print out data for time course of choice probabilities at different correlation levels TU 4/1/03
% output2 = 0;
% if (output2)
%     outfile2 = [BASE_PATH 'ProtocolSpecific\DepthDiscrim\CPtime_course_lowcorr.dat'];
%     printflag = 0;
%     if (exist(outfile2, 'file') == 0)    %file does not yet exist
%         printflag = 1;
%     end
%     fid = fopen(outfile2, 'a');
%     if (printflag)
%         fprintf(fid, 'FILE\t Nthr100\t Nthr200\t Nthr300\t Nthr400\t Nthr500\t Nthr600\t Nthr700\t Nthr800\t Nthr900\t Nthr1000\t Nthr1100\t Nthr1200\t Nthr1300\t Nthr1400\t Nthr1500\t ');
%         fprintf(fid, '\r\n');
%         printflag = 0;
%     end
% 
%     buff = sprintf('%s\t ', FILE);
%     for i=1:length(choice_prob_tc_norm)
%         buff = sprintf('%s\t %6.4f\t', buff, choice_prob_tc(i,1,1));
%     end
%     fprintf(fid, '%s', buff);
%     fprintf(fid, '\r\n');
%     buff = sprintf('%s\t ', FILE);
%     for i=1:length(choice_prob_tc_norm)
%         buff = sprintf('%s\t %6.4f\t', buff, choice_prob_tc(i,1,2));
%     end
%     fprintf(fid, '%s', buff);
%     fprintf(fid, '\r\n');
%     fclose(fid);
% end
% if (output2)
%     outfile2 = [BASE_PATH 'ProtocolSpecific\DepthDiscrim\CPtime_course_intermediatecorr.dat'];
%     printflag = 0;
%     if (exist(outfile2, 'file') == 0)    %file does not yet exist
%         printflag = 1;
%     end
%     fid = fopen(outfile2, 'a');
%     if (printflag)
%         fprintf(fid, 'FILE\t Nthr100\t Nthr200\t Nthr300\t Nthr400\t Nthr500\t Nthr600\t Nthr700\t Nthr800\t Nthr900\t Nthr1000\t Nthr1100\t Nthr1200\t Nthr1300\t Nthr1400\t Nthr1500\t ');
%         fprintf(fid, '\r\n');
%         printflag = 0;
%     end
% 
%     if(length(unique_bin_corr)>=4)
%         buff = sprintf('%s\t ', FILE);
%         for i=1:length(choice_prob_tc_norm)
%             buff = sprintf('%s\t %6.4f\t', buff, choice_prob_tc(i,4,1));
%         end
%         fprintf(fid, '%s', buff);
%         fprintf(fid, '\r\n');
%         buff = sprintf('%s\t ', FILE);
%         for i=1:length(choice_prob_tc_norm)
%             buff = sprintf('%s\t %6.4f\t', buff, choice_prob_tc(i,4,2));
%         end
%         fprintf(fid, '%s', buff);
%         fprintf(fid, '\r\n');
%     end
%     fclose(fid);
% end
% if (output2)
%     outfile2 = [BASE_PATH 'ProtocolSpecific\DepthDiscrim\CPtime_course_highcorr.dat'];
%     printflag = 0;
%     if (exist(outfile2, 'file') == 0)    %file does not yet exist
%         printflag = 1;
%     end
%     fid = fopen(outfile2, 'a');
%     if (printflag)
%         fprintf(fid, 'FILE\t Nthr100\t Nthr200\t Nthr300\t Nthr400\t Nthr500\t Nthr600\t Nthr700\t Nthr800\t Nthr900\t Nthr1000\t Nthr1100\t Nthr1200\t Nthr1300\t Nthr1400\t Nthr1500\t ');
%         fprintf(fid, '\r\n');
%         printflag = 0;
%     end
% 
%     if(length(unique_bin_corr)>=6)
%         buff = sprintf('%s\t ', FILE);
%         for i=1:length(choice_prob_tc_norm)
%             buff = sprintf('%s\t %6.4f\t', buff, choice_prob_tc(i,6,1));
%         end
%         fprintf(fid, '%s', buff);
%         fprintf(fid, '\r\n');
%         buff = sprintf('%s\t ', FILE);
%         for i=1:length(choice_prob_tc_norm)
%             buff = sprintf('%s\t %6.4f\t', buff, choice_prob_tc(i,6,2));
%         end
%         fprintf(fid, '%s', buff);
%         fprintf(fid, '\r\n');
%     end
%     fclose(fid);
% end
% 
% %-----------------------------------------------------------------------------------------------------------------------------------------------------------
% %now print out data for time course of firing for each choice separately TU 1/31/01
% output3 = 0;
% if (output3)
%     outfile3 = [BASE_PATH 'ProtocolSpecific\DepthDiscrim\PrefChoiceTime_course.dat'];
%     printflag = 0;
%     if (exist(outfile3, 'file') == 0)    %file does not yet exist
%         printflag = 1;
%     end
%     fid = fopen(outfile3, 'a');
%     if (printflag)
%         fprintf(fid, 'FILE\t Nthr100\t Nthr200\t Nthr300\t Nthr400\t Nthr500\t Nthr600\t Nthr700\t Nthr800\t Nthr900\t Nthr1000\t Nthr1100\t Nthr1200\t Nthr1300\t Nthr1400\t Nthr1500\t ');
%         fprintf(fid, '\r\n');
%         printflag = 0;
%     end
% 
%     buff = sprintf('%s\t ', FILE);
%     for i=1:length(mean_norm_pref)
%         buff = sprintf('%s\t %6.4f\t', buff, mean_norm_pref(i));
%     end
%     fprintf(fid, '%s', buff);
%     fprintf(fid, '\r\n');
%     fclose(fid);
% 
%     outfile4 = [BASE_PATH 'ProtocolSpecific\DepthDiscrim\NullChoiceTime_course.dat'];
%     printflag = 0;
%     if (exist(outfile4, 'file') == 0)    %file does not yet exist
%         printflag = 1;
%     end
%     fid = fopen(outfile4, 'a');
%     if (printflag)
%         fprintf(fid, 'FILE\t Nthr100\t Nthr200\t Nthr300\t Nthr400\t Nthr500\t Nthr600\t Nthr700\t Nthr800\t Nthr900\t Nthr1000\t Nthr1100\t Nthr1200\t Nthr1300\t Nthr1400\t Nthr1500\t ');
%         fprintf(fid, '\r\n');
%         printflag = 0;
%     end
% 
%     buff = sprintf('%s\t ', FILE);
%     for i=1:length(mean_norm_null)
%         buff = sprintf('%s\t %6.4f\t', buff, mean_norm_null(i));
%     end
%     fprintf(fid, '%s', buff);
%     fprintf(fid, '\r\n');
%     fclose(fid);
% end
% 
% %------------------------------------------------------------------------------------------------------------------
% % write out CP for each correlation level. Only use when the calculation of CP at each correlation level
% % is done with ROC_significance_test not rocN  TU 08/09/01
% output1 = 0;
% if (output1)
%     outfile1 = [BASE_PATH 'ProtocolSpecific\DepthDiscrim\CP_corr.dat'];
% 
%     printflag = 0;
%     if (exist(outfile1, 'file') == 0)    %file does not yet exist
%         printflag = 1;
%     end
%     fid = fopen(outfile1, 'a');
%     if (printflag)
%         fprintf(fid, 'File SignedBinCorr PercentPref CP CP_Pval');
%         fprintf(fid, '\r\n');
%     end
%         
%     for i = 1:length(unique_bin_corr)
%         for j = 1:length(unique_hdisp)
%             sign = (unique_hdisp(j) == Pref_HDisp)*2 - 1;	%=1 if preferred disparity, -1 if null disparity
%             signed_corr = unique_bin_corr(i) * sign;
%             percent_pref = (length(pref_dist{i,j}) / (length(pref_dist{i,j})+length(null_dist{i,j})));
% 
%             if ( (length(pref_dist{i,j}) > 0) & (length(null_dist{i,j}) > 0) )
%                 outstr1 = sprintf('%s %8.4f %8.6f %8.6f %8.6f', FILE, signed_corr, percent_pref, choice_prob(i, j), choice_prob_Pval(i, j));
%                 fprintf(fid, '%s', outstr1);
%                 fprintf(fid, '\r\n');
%             end
%         end
%     end        
%     fclose(fid);
% end
% 
% %------------------------------------------------------------------------------------------------------------------
% % write out preferred and null CP for each correlation level. Only use when the calculation of CP at each correlation level
% % is done with ROC_significance_test not rocN  TU 09/08/02
% output1 = 0;
% if (output1)
%     outfile1 = [BASE_PATH 'ProtocolSpecific\DepthDiscrim\CP_prefnull.dat'];
% 
%     printflag = 0;
%     if (exist(outfile1, 'file') == 0)    %file does not yet exist
%         printflag = 1;
%     end
%     fid = fopen(outfile1, 'a');
%     if (printflag)
%         fprintf(fid, 'File BinCorr CPpref CPnull');
%         fprintf(fid, '\r\n');
%     end
%         
%     for i = 1:length(unique_bin_corr)
%         pref_index = find(unique_hdisp == Pref_HDisp);
%         null_index = find(unique_hdisp ~= Pref_HDisp);
%         percent_pref1 = (length(pref_dist{i,pref_index}) / (length(pref_dist{i,pref_index})+length(null_dist{i,pref_index})));
%         percent_pref2 = (length(pref_dist{i,null_index}) / (length(pref_dist{i,null_index})+length(null_dist{i,null_index})));
% 
%         if ( (length(pref_dist{i,pref_index}) > 0) & (length(null_dist{i,pref_index}) > 0) & (length(pref_dist{i,null_index}) > 0) & (length(null_dist{i,null_index}) > 0))
%             outstr1 = sprintf('%s %8.4f %8.6f %8.6f %8.6f %8.6f %8.6f', FILE, unique_bin_corr(i), percent_pref1, choice_prob(i, pref_index), percent_pref2, choice_prob(i, null_index));
%             fprintf(fid, '%s', outstr1);
%             fprintf(fid, '\r\n');
%         end
%     end        
%     fclose(fid);
% end
% 
% return;
% 
% 
% %-----------------------------------------------------------------------------------------------------------------
% % write out data to analyze Z-scored spikes and vergence angle with Origin. Only use when CORRECT_FOR_VERGENCE = 1.   TU 03/14/02
% outfile = [BASE_PATH 'ProtocolSpecific\DepthDiscrim\CPverg_example.dat'];
% 
% printflag = 0;
% if (exist(outfile, 'file') == 0)    %file does not yet exist
%     printflag = 1;
% end
% fid = fopen(outfile, 'a');
% if (printflag)
%     fprintf(fid, 'Trial VergAngle ZSpike Residual Choice ');
%     fprintf(fid, '\r\n');
% end
%         
% for i = 1:length(binoc_corr)
%             outstr1 = sprintf('%8.6f %8.6f %8.6f %8.6f %8.6f', trials(i), calib_h_verg(i), Z_Spikes(i), r(i), choice(i));
%             fprintf(fid, '%s', outstr1);
%             fprintf(fid, '\r\n');
% end
% 
% fclose(fid);
% 
% %-----------------------------------------------------------------------------------------------------------------
% % write out data to analyze Z-scored spikes with Origin. Only use when CORRECT_FOR_SLOW_SPIKE_RATE_CHANGE = 1.   TU 08/09/01
% i = size(PATH,2) - 1;
% while PATH(i) ~='\'	%Analysis directory is one branch below Raw Data Dir
%     i = i - 1;
% end   
% PATHOUT = [PATH(1:i) 'Analysis\Temp\'];
% outfile1 = [PATHOUT 'ZSpike.dat'];
% 
% printflag = 0;
% if (exist(outfile1, 'file') == 0)    %file does not yet exist
%     printflag = 1;
% end
% fid = fopen(outfile1, 'a');
% if (printflag)
%     fprintf(fid, 'ZSpike Trial SignedBinCorr Choice');
%     fprintf(fid, '\r\n');
% end
%         
% for i = 1:length(binoc_corr)
%             outstr1 = sprintf('%8.4f %8.6f %8.6f %8.6f', Z_Spikes(i), trials(i), signed_bin_corr(i), choice(i));
%             fprintf(fid, '%s', outstr1);
%             fprintf(fid, '\r\n');
% end
% 
% fclose(fid);
% 
% %-----------------------------------------------------------------------------------------------------------------
% % write out data to analyze Z-scored spikes with Origin. Only use when CORRECT_FOR_SLOW_SPIKE_RATE_CHANGE = 1.   TU 08/09/01
% outfile1 = [BASE_PATH 'ProtocolSpecific\DepthDiscrim\ZSpike.dat'];
% 
% printflag = 0;
% if (exist(outfile1, 'file') == 0)    %file does not yet exist
%     printflag = 1;
% end
% fid = fopen(outfile1, 'a');
% if (printflag)
%     fprintf(fid, 'B BINT BINT PValue');
%     fprintf(fid, '\r\n');
% end
%         
% outstr1 = sprintf('%8.6f %8.6f %8.6f %8.6f', b(1), bint(1,1), bint(1,2), stats(3));
% fprintf(fid, '%s', outstr1);
% fprintf(fid, '\r\n');
% 
% fclose(fid);