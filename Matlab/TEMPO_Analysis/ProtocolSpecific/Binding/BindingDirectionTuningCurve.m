%-----------------------------------------------------------------------------------------------------------------------
%-- BindingDirectionTuningCurve.m -- Plots Direction tuning curve and computes/plots a Gaussian fits to these curves
%--    Works on both spike trains and LFPs.
%--	Adapted from DirectionTuningCurve.m by BJP
%-----------------------------------------------------------------------------------------------------------------------
function BindingDirectionTuningCurve(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE);

% not implemented yet to select output
output = 1;

TEMPO_Defs; %contains protocol specific keywords - 1/4/01 BJP

switch SpikeChan
    case 4 %first LFP channel
        neural_db = LFP_DB;
        SpikeChan = 1;
    case 5 %second LFP channel
        neural_db = LFP_DB;
        SpikeChan = 2;
    case 1 %first Spike Channel
        neural_db = SPIKE_DB;
        SpikeChan = 1;
    case 2 %second Spike Channel
        neural_db = SPIKE_DB;
        SpikeChan = 2;
end



symbols = {'bo' 'r*' 'go' 'mo' 'b*' 'r*' 'g*' 'c*'};
lines = {'b-' 'r-.' 'g-' 'm-' 'b--' 'r--' 'g--' 'c--'};

curr_fig = figure;
set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [150 100 500 573], 'Name', 'Direction Tuning Curve');
subplot(2, 1, 2);


direction = data.obj_params(OBJ_TRAJ_ORIENT,:,1);
       
%get indices of any NULL conditions (for measuring spontaneous activity
%null_trials = logical( (direction == data.one_time_params(NULL_VALUE)) );
null_trials = logical( (data.misc_params(CONDITION,:) == 0) );  		%null condition = condition 0

%now, remove trials from direction and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:length(direction);		% a vector of trial indices
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );

%   [direction' null_trials' select_trials']

N_peaks = 1; %2 for V1 cells

if (Protocol == BIND_DIR_TUNING)
    dir(1,:) = data.obj_params(OBJ_TRAJ_ORIENT,:,1) == 1;
    dir(2,:) = data.obj_params(OBJ_TRAJ_ORIENT,:,1) == -1;

    if (N_peaks == 1)
        direction(dir(2,:)) = 180 + direction(dir(2,:));
        dir(1,:) = dir(1,:) + dir(2,:);
    end
else 
    dir(1,:) = data.obj_params(OBJ_DIR,:,1) == 1;
    
    auto_channel_select = 1;
    

    switch SpikeChan
    case 1
        select_trials = select_trials & (data.obj_params(OBJ_STATUS,:,1) == 1) | null_trials;
    case 2
        select_trials = select_trials & (data.obj_params(OBJ_STATUS,:,1) == 0) | null_trials;
    end
    if sum(select_trials == 0)
        %find which channel we are mapping and override current channel
        %being analyzed.  Leave alone if channel correct or mapping both
        select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );
        if auto_channel_select == 1

            SpikeChan = ~(SpikeChan - 1) + 1;   % switch spike channels
        end    
        select_trials = select_trials & (data.obj_params(OBJ_STATUS,:,SpikeChan) == 1) | null_trials;
    end   
end
    
switch neural_db
    case SPIKE_DB
        %now, get the firing rates for all the trials 
        spike_rates = data.spike_rates(SpikeChan, :);
    case LFP_DB
       LFP_bin_width = (data.htb_header{LFP_DB}.skip + 1) / (data.htb_header{LFP_DB}.speed_units / data.htb_header{LFP_DB}.speed);
       for trial = 1:size(data.lfp_data,3);
            lfp = data.lfp_data(SpikeChan, floor(StartEventBin(trial)/(LFP_bin_width*1000)) :floor(StopEventBin( trial )/(LFP_bin_width*1000) ), trial);
            [freq, ampl] = FourierTransform_1D(linspace(0,LFP_bin_width*length(lfp),length(lfp)), lfp, length(lfp), 1,0);
            %use the power of frequencies between 100 and 200 hz
            spike_rates(trial) = sum( (ampl( (freq > 100) & (freq <= 200) )  ).^2 ); 
        end
    end        
        

for set = 1:N_peaks

	plot_x = direction(~null_trials & select_trials & dir(set, :) );
	plot_y = spike_rates(~null_trials & select_trials &  dir(set, :) );


	%NOTE: inputs to PlotTuningCurve must be column vectors, not row vectors, because of use of munique()
	[px, py, perr] = PlotTuningCurve(plot_x', plot_y', symbols{set}, lines{set}, 1, 0);	%note: last arg=0 means just get output, no plot

	unique_dirs = px; % store category groups for ANOVAs

	%now, shift the px, py, and perr vectors such that the peak of tuning curve is in middle of axis range
	% now, we need to shift these vectors so that the peak response is always in middle of vector
	ctr_indx = round(length(px)/2 - rem(length(px),2)) + 1;
	[max_val max_indx] = max(py);
   shift = max_indx - ctr_indx;

   if (N_peaks == 1)
    if (shift > 0)
    	px = [ px(shift+1 : length(px)) ; px(1 : shift)+360];
     	py = [ py(shift+1 : length(py)) ; py(1 : shift)];
       perr = [ perr(shift+1 : length(perr)) ; perr(1 : shift)];
   end
 	if (shift < 0)
    	px = [ px(length(px)+shift+1 : length(px))-360 ; px(1 : length(px)+shift)];
     	py = [ py(length(py)+shift+1 : length(py)) ;  py(1 : length(py)+shift)];
     	perr = [ perr(length(perr)+shift+1 : length(perr)) ;  perr(1 : length(perr)+shift)];
    end


	% Since direction is circular, duplicate  the lowest value (px) as px+ 360, and change spike arrays accordingly
	 px = [px; px(1)+360];
	 py = [py; py(1)];
	 perr = [perr; perr(1)];
 end
	hold on;
	%plot the data, after shifting as necessary above
	errorbar(px, py, perr, perr, symbols{set} );
	hold on;

	%now, fit the data with a Gaussian curve and plot this as well
	means = [px py];
	raw = [px py];

	[pars] = gaussfit(means, raw, 1);
	x_interp = (px(1)): 1 : (px(length(px)));
	y_interp = gaussfunc(x_interp, pars);
   yl(set) = max(y_interp');
   ym(set) = min(y_interp');
   % Note: this func MUST match that in gaussfit.m
	h(set) = plot(x_interp, y_interp, lines{set});

	%now, get the firing rate for NULL condition trials and add spontaneous rate to plot
	null_x = [min(px) max(px)];
	null_rate = mean(spike_rates(null_trials & select_trials) );
	null_y = [null_rate null_rate];
	hold on;
	plot(null_x, null_y, 'k--');
	hold off;



	% calculate some metrics and stats then print them in plot
	pref_dir(set) = pars(3);
	p_value(set) = spk_anova(plot_y, plot_x, unique_dirs);
	base_rate(set) = pars(1);
	amplitude(set) = pars(2);
	max_rate(set) = base_rate(set) + amplitude(set);
	width(set) = sqrt(-(log(.5)))*pars(4)*2*sqrt(2);
	DSI(set) = 1 - (base_rate(set) - null_rate)/(max_rate(set) - null_rate); 

end

yl = max(yl);
%yl = max(raw(:,2));
ym = min(raw(:,2));

YLim([0 yl] );
        
XLabel('Direction of Motion (deg)');
YLabel('Response (spikes/sec)');

% Compute a direction discrimination index analogous to the DDI
[DirDI, var_term] = Compute_DDI(plot_x, plot_y);


%now, print out some useful information in the upper subplot
subplot(2, 1, 1);
PrintGeneralData(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

axis([0 100 0 100]);
axis('off');
font_size = 9;
bump_size = 7;

xpos = -10;   
ypos = 20;

line = sprintf('Fitted Tuning Curve:');
text(xpos,ypos,line,'FontSize',font_size+2);		ypos = ypos - bump_size;      
line = sprintf('Preferred Direction = %1.5g   Resp = %0.5g', pref_dir, max_rate);
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;    
line = sprintf('Base Rate = %0.5g   Spont = %0.5g', base_rate, null_rate);
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('Amplitude = %0.3g', amplitude);
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;       
line = sprintf('Width (FWHM) = %0.3g', width);
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;         
line = sprintf('DSI = %0.3g', DSI);
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;       

xpos = 60;   
ypos = 10;
font_size = 8;
line = sprintf('ANOVA: P = %0.3g', p_value);
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;       
line = sprintf('Dir Discrim Ind = %0.5g', DirDI);
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;  

if N_peaks > 1
    legend(h,'0<dir<180','180<dir<360', 4);
end


    %output tuning curve metrics
    if (output == 1)
        i = size(PATH,2) - 1;
        while PATH(i) ~='\'	%Analysis directory is one branch below Raw Data Dir
            i = i - 1;
        end   
        PATHOUT = [PATH(1:i) 'Analysis\Tuning\'];
        i = size(FILE,2) - 1;
        while FILE(i) ~='.'
            i = i - 1;
        end
        FILEOUT = [FILE(1:i) 'dir'];
        eval(['save ' PATHOUT FILEOUT  ' pref_dir max_rate base_rate null_rate amplitude width DSI p_value DirDI SpikeChan StartCode StopCode BegTrial EndTrial StartOffset StopOffset StartEventBin StopEventBin PATH FILE'])
        
        %save figure
        FILEOUT = [FILE(1:i) 'fig'];        
        saveas (curr_fig, [PATHOUT FILEOUT]);
    
    end
return;