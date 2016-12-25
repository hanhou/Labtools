%-----------------------------------------------------------------------------------------------------------------------
%-- CombineDirecDiscrim.m -- Does not plot anything but instead saves the
%performance data to a .mat file where the data can be accumulated.  Then
%the data is combined across all data (across multiple runs) in that .mat file and is used to
%compute the overall percentages and threshold (across multiple runs).
%Critically, there must exist a file DirectionDiscrim\CumulativeFiles\DirecDataSort.txt
%which contains lines formmated 'FILE injection day', where FILE is the
%filename, injection is the injection #, and day is the day number:
%1=pre, 2=1hr post, 3=1d post, 4=2d post, 5=3d post
%--	VR, 6/15/06
%-----------------------------------------------------------------------------------------------------------------------
function [monkey_alpha] = CombineDirecDiscrim(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

TEMPO_Defs;		
Path_Defs;
ProtocolDefs;	%needed for all protocol specific functions - contains keywords - BJP 1/4/01

%get the column of values of directions in the dots_params matrix
direction = data.dots_params(DOTS_DIREC,:,PATCH1);
unique_direction = munique(direction');

%get the patch X-center location
x_ctr = data.dots_params(DOTS_AP_XCTR, :, PATCH1);
unique_x_ctr = munique(x_ctr');

%get the motion coherences
coherence = data.dots_params(DOTS_COHER, :, PATCH1);
unique_coherence = munique(coherence');

%now, get the firing rates for all the trials 
spike_rates = data.spike_rates(SpikeChan, :);

%get indices of any NULL conditions (for measuring spontaneous activity
null_trials = logical( (coherence == data.one_time_params(NULL_VALUE)) );

%now, select trials that fall between BegTrial and EndTrial
trials = 1:length(coherence);
%a vector of trial indices
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );

[direction' coherence' spike_rates' null_trials' select_trials'];

Pref_direction = data.one_time_params(PREFERRED_DIRECTION);



%% *********** PSYCHOMETRIC ANALYSIS ****************************

pct_correct = []; N_obs = []; fit_data = [];
monkey_alpha = []; monkey_beta = [];

%first compute correct_trials and total_trials
for j = 1:length(unique_x_ctr)
    for i=1:length(unique_coherence)
        trial_inds = ((x_ctr == unique_x_ctr(j)) & (coherence == unique_coherence(i)) & select_trials);
        correct_trials_inds = (trial_inds & (data.misc_params(OUTCOME, :) == CORRECT) );
        corr_trials(j,i) = sum(correct_trials_inds);
        total_trials(j,i) = sum(trial_inds);
        pct_correct(j,i) = corr_trials(j,i)/total_trials(j,i);
    end
end

%next, extract the condition params from the file that Syed made.
infile = [BASE_PATH 'ProtocolSpecific\DirectionDiscrim\CumulativeFiles\DirecDataSort.txt'];
[AllFiles,AllInjections,AllDays] = textread(infile,'%s%d%d');
file_ind = find(strcmp(AllFiles,FILE),1);
if (length(file_ind) == 0)    
    disp('File Not Found in DirecDataSort.txt');
    return
end
inj = AllInjections(file_ind);
day = AllDays(file_ind);

%then, append these data to the appropriately named file.
outfile = sprintf('%sProtocolSpecific\\DirectionDiscrim\\CumulativeFiles\\Direc_Inj%d_Day%d.mat',BASE_PATH,inj,day);
if (exist(outfile, 'file') == 0)
    %output = {unique_x_ctr unique_coherence corr_trials total_trials};
    cum_corr_trials = corr_trials;
    cum_total_trials = total_trials;
    save(outfile, 'unique_x_ctr', 'unique_coherence', 'cum_corr_trials', 'cum_total_trials');
else
    load(outfile);
    cum_corr_trials = cum_corr_trials + corr_trials;
    cum_total_trials = cum_total_trials + total_trials;
    pct_correct = cum_corr_trials./cum_total_trials;     %need to update pct_correct for the rest of the protocol
    save(outfile, 'unique_x_ctr', 'unique_coherence', 'cum_corr_trials', 'cum_total_trials');
end

%now continue and compute threshold
for j = 1:length(unique_x_ctr)
    for i=1:length(unique_coherence)
        % data for Weibull fit
        fit_data(j,i,1) = unique_coherence(i);
        fit_data(j,i,2) = pct_correct(j,i);
        fit_data(j,i,3) = cum_total_trials(j,i);
    end
    [monkey_alpha(j) monkey_beta(j)]= weibull_fit(squeeze(fit_data(j,:,:)));
end

%     line = sprintf('Monkey: Xctr = %6.2f threshold = %6.3f %%, slope = %6.3f', unique_x_ctr(j), monkey_alpha(j), monkey_beta(j) );
%     pct_correct(1), pct_correct(2), pct_correct(3), pct_correct(4), pct_correct(5), pct_correct(6) );

output = 1;
if (output)
    
    %------------------------------------------------------------------------
    %write out all relevant parameters to a cumulative text file, VR 6/15/06
    outfile = [BASE_PATH 'ProtocolSpecific\DirectionDiscrim\Combined_psycho_Curve_summary.dat'];
    printflag = 0;
    if (exist(outfile, 'file') == 0)    %file does not yet exist
        printflag = 1;
    end
    fid = fopen(outfile, 'a');
    if (printflag)
        fprintf(fid, 'FILE\t Inj\t Day\t PrDir\t PrSpd\t PrHDsp\t RFX\t RFY\t RFDiam\t MThr\t MSlp\t Ntrials\t Coher1\t Pct1\t Coher2\t Pct2\t Coher3\t Pct3\t Coher4\t Pct4\t Coher5\t Pct5\t Coher6\t Pct6\t ');
        fprintf(fid, '\r\n');
        printflag = 0;
    end

    for j = 1:length(unique_x_ctr)
        buff = sprintf('%s\t %d\t %d\t %6.1f\t %6.2f\t %6.3f\t %6.2f\t %6.2f\t %6.2f\t %6.3f\t %6.4f\t %6.3f\t ', ...
            FILE, inj, day, data.neuron_params(PREFERRED_DIRECTION, 1), data.neuron_params(PREFERRED_SPEED, 1), data.neuron_params(PREFERRED_HDISP, 1), ...
            unique_x_ctr(j), data.neuron_params(RF_YCTR, 1), data.neuron_params(RF_DIAMETER, 1),...
            monkey_alpha(j),monkey_beta(j), sum(cum_total_trials(j,:)) );
        for i = 1:length(unique_coherence)
            buff = [buff sprintf('%6.2f\t %6.3f\t ',unique_coherence(i),pct_correct(j,i))];
        end
        fprintf(fid, '%s', buff);
        fprintf(fid, '\r\n');
    end
    
    fclose(fid);
    %------------------------------------------------------------------------
end

return;