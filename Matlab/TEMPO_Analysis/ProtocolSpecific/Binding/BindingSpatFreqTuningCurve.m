%-----------------------------------------------------------------------------------------------------------------------
%-- DirectionTuningCurve2.m -- Plots 2 direction tuning curves and computes/plots a Gaussian fits to these curves
%--    
%--	Adapted from DirectionTuningCurve.m by BJP
%-----------------------------------------------------------------------------------------------------------------------
function BindingSpatFreqTuningCurve(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

% not implemented yet to select output
output = 1;

ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP

symbols = {'bo' 'r*' 'go' 'mo' 'b*' 'r*' 'g*' 'c*'};
lines = {'b-' 'r-.' 'g-' 'm-' 'b--' 'r--' 'g--' 'c--'};

spat_freq = data.obj_params(OBJ_AMPL,:,1);

unique_freq = munique(spat_freq);
%now, get the firing rates for all the trials 
spike_rates = data.spike_rates(SpikeChan, :);

    %some stuff done to send data to Izumi: GCD, 2/21/01
    %disp('saving...');
    %spikes_mat = squeeze(data.spike_data(1,:,:))';
    %size(spikes_mat)
    %outfid = fopen('temp2.dat', 'w');
    %for i=1:size(spikes_mat, 1)
    %    fprintf(outfid, '%d ', find(spikes_mat(i,:)>0) );
    %    fprintf(outfid, '\n');
    %end
    %fclose(outfid);


%get indices of any NULL conditions (for measuring spontaneous activity
null_trials = logical( (data.misc_params(CONDITION,:) == 0) );  		%null condition = condition 0
%null_trials = logical( (direction == data.one_time_params(NULL_VALUE)) );

%now, remove trials from direction and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:length(spat_freq);		% a vector of trial indices
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );

%   [spat_freq' null_trials' select_trials']

figure;
set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [150 100 500 573], 'Name', 'spat_freq Tuning Curve');
subplot(2, 1, 2);


	plot_x = spat_freq(~null_trials & select_trials );
	plot_y = spike_rates(~null_trials & select_trials );


	%NOTE: inputs to PlotTuningCurve must be column vectors, not row vectors, because of use of munique()
	[px, py, perr] = PlotTuningCurve(plot_x', plot_y', symbols, lines, 1, 0);	%note: last arg=0 means just get output, no plot

 
	% Since spat_freq is circular, duplicate  the lowest value (px) as px+ 360, and change spike arrays accordingly
	% px = [px; px(1)+360];
	% py = [py; py(1)];
	% perr = [perr; perr(1)];

	hold on;
	%plot the data, after shifting as necessary above
	errorbar(px, py, perr, perr, symbols{1} );
	hold on;

	%now, fit the data with a Gaussian curve and plot this as well
	means = [px py];
	raw = [px py];

	[pars] = gammafit(means, raw);
	x_interp = (px(1)): 0.1 : (px(length(px)));
	y_interp = gammafunc(x_interp, pars);
    yl = max(y_interp');
   
	h = plot(x_interp, real(y_interp), lines{1});

	%now, get the firing rate for NULL condition trials and add spontaneous rate to plot
	null_x = [min(px) max(px)];
	null_rate = mean(data.spike_rates(SpikeChan, (null_trials & select_trials)  ));
	null_y = [null_rate null_rate];
	hold on;
	plot(null_x, null_y, 'k--');
	hold off;



	% calculate some metrics and stats then print them in plot
	pref_freq = pars(3);
	p_value = spk_anova(plot_y, plot_x, unique_freq);
	base_rate = pars(1);
	amplitude = pars(2);
	max_rate = base_rate + amplitude;
	width = sqrt(-(log(.5)))*pars(4)*2*sqrt(2);

yl = max(yl);
YLim([0 yl]);	% set the lower limit of the Y axis to zero
Xlim([0 max(x_interp)]);
XLabel('Spatial Frequency (c/deg)');
YLabel('Response (spikes/sec)');


%now, print out some useful information in the upper subplot
subplot(2, 1, 1);
PrintGeneralData(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

PrintFrequencyData(p_value, base_rate, null_rate, amplitude, pref_freq, max_rate, width); 

return;