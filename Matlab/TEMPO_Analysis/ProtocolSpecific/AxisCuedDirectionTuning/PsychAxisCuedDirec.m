%-----------------------------------------------------------------------------------------------------------------------
%-- PsychAxisCuedDirec.m -- Plots direction tuning curve with and without cues
%--	VR, 5/24/05
%-----------------------------------------------------------------------------------------------------------------------
function PsychAxisCuedDirec(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE);

TEMPO_Defs;		
Path_Defs;
ProtocolDefs;	%needed for all protocol specific functions - contains keywords - BJP 1/4/01

%get the column of values of target directions in the dots_params matrix
targ_direc = data.cue_params(AXIS_CUE_DIREC,:,PATCH2);
unique_targ_direc = munique(targ_direc');

%get the motion coherences
coherence = data.dots_params(DOTS_COHER, :, PATCH1);
unique_coherence = munique(coherence');

%get the cue status for each trial.
cue_status = data.cue_params(AXIS_CUE_STATUS, :, PATCH2);
unique_cue_status = munique(cue_status');

%get the deviation of motion directions around the targets for each trial.
direc_spread = data.cue_params(DIREC_SPREAD, :, PATCH1);
unique_direc_spread = munique(direc_spread');

%now, get the firing rates for all the trials 
spike_rates = data.spike_rates(SpikeChan, :);

%get indices of any NULL conditions (for measuring spontaneous activity)
null_trials = logical( (coherence == data.one_time_params(NULL_VALUE)) );

%now, select trials that fall between BegTrial and EndTrial
trials = 1:length(coherence);
%a vector of trial indices
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );

%[direction' coherence' spike_rates' null_trials' select_trials']

Pref_direction = data.one_time_params(PREFERRED_DIRECTION);
figure;
set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [250 50 500 573], 'Name', 'Performance vs Direc Spread');
%subplot(2, 1, 2);

symbols = {'ko', 'kx', 'bo', 'rx'};
lines = {'k-', 'k--', 'b:', 'r:'};

%% ************* ANALYZE PERFORMANCE **************************
%for each coherence, for each cue status, for each spread, compute percent
%correct and plot
pct_correct = []; N_obs = []; 

for i = 1:length(unique_coherence)
    for j = 1:length(unique_targ_direc)/2
        pct_cue =[];
        for k = 1:length(unique_cue_status)
            for l = 1:length(unique_direc_spread)
                rel_trials = logical ( (coherence == unique_coherence(i)) & ...
                    ((targ_direc == unique_targ_direc(j)) | (targ_direc == unique_targ_direc(j+length(unique_targ_direc)/2)) ) & ...
                    (cue_status == unique_cue_status(k)) & ...
                    (direc_spread == unique_direc_spread(l)) );
                %rel trials includes trials going in both directions
                correct_trials = (rel_trials & (data.misc_params(OUTCOME,:) == CORRECT));
                pct_correct(i,j,k,l) = sum(correct_trials)/sum(rel_trials);
                N_obs(i,j,k,l) = sum(rel_trials);
            end
            pct_cue(:,k) = pct_correct(i,j,k,:);
        end
        subplot(length(unique_coherence),length(unique_targ_direc)/2, (i-1)*length(unique_targ_direc)/2 + j);
        bar(unique_direc_spread, pct_cue);
        colormap([0 0 1; 1 1 0]); %set colors to blue and yellow for easy b/w printing...
        xlabel('Direc Spread');
        YLim([0,1])
        if(j == 1) %first column
            ylabel(sprintf('Coher = %d%%\nPct Correct',unique_coherence(i)));
        else
            ylabel('Pct Correct');
        end
        
        if (i == 1) %top row
            title(sprintf('%d/%d axis',unique_targ_direc(j),unique_targ_direc(j+length(unique_targ_direc)/2)));
        end
    end
end

% subplot(length(unique_coherence),length(unique_targ_direc)/2,1);
% if (length(unique_cue_status) == 2)
%     legend('Cue Off','Cue On');
% else
%     if(unique_cue_status(k) == 1)
%         legend ('Cue On');
%     else
%         legend ('Cue Off');
%     end
% end

% now compute performance histographs collapsing across coherence and the
% two axes.
figure;
pct_cue = []; pct_correct = [];
for k = 1:length(unique_cue_status)
    for l = 1:length(unique_direc_spread)
        rel_trials = logical ( (cue_status == unique_cue_status(k)) & ...
            (direc_spread == unique_direc_spread(l)) );
        %rel trials includes trials going in both directions
        correct_trials = (rel_trials & (data.misc_params(OUTCOME,:) == CORRECT));
        pct_correct(k,l) = sum(correct_trials)/sum(rel_trials);
        N_obs(k,l) = sum(rel_trials);
    end
    pct_cue(:,k) = pct_correct(k,:);
end
bar(unique_direc_spread, pct_cue);
colormap([0 0 1; 1 1 0]);
xlabel('Direc Spread');
ylabel('Pct Correct');
YLim([0,1]);

%% *********** PSYCHOMETRIC ANALYSIS ****************************
figure;
set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [250 50 500 573], 'Name', 'Psychometric Function');
subplot(2, 1, 2);

pct_correct = []; N_obs = []; fit_data = [];
monkey_alpha = []; monkey_beta = [];

%computes psychometric curves for both axis directions and cue status
for j = 1:length(unique_targ_direc)/2 + length(unique_cue_status)
    for i=1:length(unique_coherence)
        if (j <= length(unique_targ_direc)/2)
            trials = ((coherence == unique_coherence(i)) & ( (targ_direc == unique_targ_direc(j)) | ...
            (targ_direc == unique_targ_direc(j + length(unique_targ_direc)/2)) ) ...
            & select_trials);
        else
            trials = ( (coherence == unique_coherence(i)) & ...
                (cue_status == unique_cue_status(j-length(unique_targ_direc)/2)) & select_trials);
        end
        correct_trials = (trials & (data.misc_params(OUTCOME, :) == CORRECT) );
        pct_correct(i) = sum(correct_trials)/sum(trials);
        N_obs(i) = sum(trials);
        % data for Weibull fit
        fit_data(i, 1) = unique_coherence(i);
        fit_data(i, 2) = pct_correct(i);
        fit_data(i, 3) = N_obs(i);
    end
    hold on;
    Handl(1) = plot(unique_coherence, pct_correct, symbols{j});
    hold off;
    
    fit_x = unique_coherence(1):0.1: unique_coherence(length(unique_coherence));
    [monkey_alpha(j) monkey_beta(j)]= weibull_fit(fit_data);
    monkey_fit_y = weibull_curve(fit_x, [monkey_alpha(j) monkey_beta(j)]);
 
    hold on;
    plot(fit_x, monkey_fit_y, lines{j});
    hold off;
end

xlabel('Coherence (% dots)');
ylabel('Fraction Correct');

YLim([0.4 1]);
%comment out the next 2 lines if you want the plot to be on a LINEAR X-axis
set(gca, 'XScale', 'log');
XLim([1 100]);

%now, print out some useful information in the upper subplot
subplot(2, 1, 1);
PrintGeneralData(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, PATH, FILE);

start_time = find(data.event_data(1, :, 1) == VSTIM_ON_CD);
stop_time = find(data.event_data(1, :, 1) == VSTIM_OFF_CD);
stim_duration = stop_time - start_time

%now, print out some specific useful info.
xpos = 0; ypos = 10;
font_size = 11;
bump_size = 8;
for j = 1:length(unique_targ_direc)/2
    line = sprintf('Monkey: Axis = %d/%d, threshold = %6.3f %%, slope = %6.3f', ...
        unique_targ_direc(j),unique_targ_direc(j+length(unique_targ_direc)/2), monkey_alpha(j), monkey_beta(j) );
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
end
for j = 1:length(unique_cue_status)
    line = sprintf('Monkey: Cue Status = %d, threshold = %6.3f %%, slope = %6.3f', ...
        unique_cue_status(j), monkey_alpha(j+length(unique_targ_direc)/2), monkey_beta(j+length(unique_targ_direc)/2) );
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
end 
% line = sprintf('Directions tested: %6.3f, %6.3f deg', unique_direction(1), unique_direction(2) );
% text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
% line = sprintf('Stimulus Duration: %5d', stim_duration );
% text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
% line = sprintf('(6,%0.5g) (12,%0.5g) (24,%0.5g) (48,%0.5g)', ...
%     pct_correct(1), pct_correct(2), pct_correct(3), pct_correct(4), pct_correct(5), pct_correct(6) );
% text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;

% output = 1;
% if (output)
%     
%     %------------------------------------------------------------------------
%     %write out all relevant parameters to a cumulative text file, GCD 8/08/01
%     outfile = [BASE_PATH 'ProtocolSpecific\DirecDiscrim\Psycho_Curve_summary.dat'];
%     printflag = 0;
%     if (exist(outfile, 'file') == 0)    %file does not yet exist
%         printflag = 1;
%     end
%     fid = fopen(outfile, 'a');
%     if (printflag)
%         fprintf(fid, 'FILE\t\t PrDir\t PrSpd\t PrHDsp\t RFX\t RFY\t RFDiam\t MThr\t MSlp\t DspLo\t DspHi\t Ntrials\t HCorr\t Durat\t');
%         fprintf(fid, '\r\n');
%         printflag = 0;
%     end
%     buff = sprintf('%s\t %6.1f\t %6.2f\t %6.3f\t %6.2f\t %6.2f\t %6.2f\t %6.3f\t %6.4f\t %6.3f\t %6.3f\t %4d\t %6.3f\t %5d\t', ...
%         FILE, data.neuron_params(PREFERRED_DIRECTION, 1), data.neuron_params(PREFERRED_SPEED, 1), data.neuron_params(PREFERRED_HDISP, 1), data.neuron_params(RF_XCTR, 1), data.neuron_params(RF_YCTR, 1), data.neuron_params(RF_DIAMETER, 1),...
%         monkey_alpha,monkey_beta,unique_direction(1), unique_direction(2), (1+ EndTrial - BegTrial), unique_coherence(length(unique_coherence)), stim_duration );
%     fprintf(fid, '%s', buff);
%     fprintf(fid, '\r\n');
%     fclose(fid);
%     %------------------------------------------------------------------------
% end

return;