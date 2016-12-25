function [perr, px, py, Uncorr_resp, plot_x, plot_y] = HDispTuningParams(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, UseSyncPulses);
%for h. disp i need to obtain the following info:
   %perr
   %px
   %py
   %uncorr_resp
   %plot_x
   %plot_y

ProtocolDefs;

%following function checks offset times and outputs starting and ending analysis (in spike bins)
num_trials = size(data.event_data, 3);
[StartOffsetBin StopOffsetBin] = CheckTimeOffset(data, num_trials, StartCode, StopCode, StartOffset, StopOffset, UseSyncPulses);

if (~isempty(data.spike_data))
   %compute the firing rate over all trials during the period between StartCode and StopCode
   data.spike_rates = ComputeSpikeRates(data, num_trials, StartCode, StopCode, StartOffsetBin, StopOffsetBin);
end

symbols = {'bo' 'ro' 'go' 'ko' 'b*' 'r*' 'g*' 'k*'};
lines = {'b-' 'r-' 'g-' 'k-' 'b--' 'r--' 'g--' 'k--'};

%get the column of values of horiz. disparity in the dots_params matrix
hor_disp = data.dots_params(DOTS_HDISP,:,PATCH1);

%get indices of any NULL conditions (for measuring spontaneous activity)
null_trials = logical( (hor_disp == data.one_time_params(NULL_VALUE)) );

%get the column of speed values
speed = data.dots_params(DOTS_SPEED,:,PATCH1);
unique_speed = munique(speed(~null_trials)');

%now, get the firing rates for all the trials 
spike_rates = data.spike_rates(SpikeChan, :);

%get indices of monoc. and uncorrelated controls
control_trials = logical( (hor_disp == LEYE_CONTROL) | (hor_disp == REYE_CONTROL) | (hor_disp == UNCORR_CONTROL) );

%now, remove trials from hor_disp and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:length(hor_disp);		% a vector of trial indices
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );

% Calculate spontaneous rates before looping through so can calculate DTI
null_rate = mean(data.spike_rates(SpikeChan, null_trials & select_trials));

for i=1:length(unique_speed)	%for each different speed value, plot a separate disparity tuning curve
    speed_select = logical( (speed == unique_speed(i)) );
    
    plot_x = hor_disp(speed_select & ~null_trials & ~control_trials & select_trials);
    plot_y = spike_rates(speed_select & ~null_trials & ~control_trials & select_trials); 

    %NOTE: inputs to PlotTuningCurve must be column vectors, not row vectors, because of use of munique()
    [px, py, perr, pmax(i), pmin(i)] = PlotTuningCurve(plot_x', plot_y', symbols{i}, lines{i}, 1, 0);
    
    Uncorr_trials = logical( (hor_disp == UNCORR_CONTROL) );
    Uncorr_resp = spike_rates(speed_select & select_trials & Uncorr_trials);
    Uncorr_resp = mean(Uncorr_resp);

end

