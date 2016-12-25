%*******************************************************************************************************************
% Coherence.m - Generates Cross Correlograms Smoothed Shuffled Correlograms, and subtracted 
%       correlograms.  Then does bootstrapping.  In the process, subtracts shuffled from raw, and subtracts mean of 
%       subtracted correlogram to set baseline of 0.  This is done for each correlation going into the shuffles and 
%       cross correlogram - analyzes in frequency domain
%       BJP 1/10/01

function Compute_Coherence(data, Protocol, Analysis, SpikeChan, SpikeChan2, start_code, stop_code, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE);

line_types = ['b-'; 'r-'; 'g-'; 'k-'; 'g-'; 'g-'; 'g-'; 'g-'];
TEMPO_Defs;
ProtocolDefs;

segment_length = 512;
step_increment = 100;
segment_start_times = 1:step_increment:3000;
segment_end_times = segment_start_times + segment_length - 1; 

spec_time = 0:step_increment:3000;

%load msac times and regions with evoked response...
plot_PSTH = 1;
if (plot_PSTH == 1)
    i = size(PATH,2) - 1;
    while PATH(i) ~='\'	%Analysis directory is one branch below Raw Data Dir
        i = i - 1;
    end   
    PATHOUT = [PATH(1:i) 'Analysis\PSTH\'];
    i = size(FILE,2) - 1;
    while FILE(i) ~='.'
        i = i - 1;
    end
    FILEOUT = [FILE(1:i) 'psh'];
    eval(['load ' PATHOUT FILEOUT ' time hist_data2 hist_data2B minimum_segment_length mini_segment_start_times mini_segment_end_times -MAT'])
    hist_time = time;
    time = [];
else
    mini_segment_start_times{1} = segment_start_times
    mini_segment_end_times{1} = segment_end_times    
    %mini_segment_end_times{cond} = mini_segment_end_times{1};
    %mini_segment_start_times{cond} = mini_segment_start_times{1};
end



%keep track of spike channels before reindexing
SpikeChanA = SpikeChan;
SpikeChanB = SpikeChan2;

num_spike_channels = size(data.spike_data,1);

switch SpikeChan
    case num_spike_channels + 1%first LFP channel
        neural_db1 = LFP_DB;
        SpikeChan = 1;
    case num_spike_channels + 2%second LFP channel
        neural_db1 = LFP_DB;
        SpikeChan = 2;
    otherwise 
        neural_db1 = SPIKE_DB;
end

%assign the proper dbase index for the neural signals and the correct indices for channel 2
%allow for the possibility of using the same channel twice -- autocorrelogram

switch SpikeChan2
    case num_spike_channels + 1%first LFP channel
        neural_db2 = LFP_DB;
        SpikeChan2 = 1;
    case num_spike_channels + 2%second LFP channel
        neural_db2 = LFP_DB;
        SpikeChan2 = 2;
    otherwise 
        neural_db2 = SPIKE_DB;
end

%first, we'll need the bin_width (in sec) of our spike raster + events log
h = data.htb_header{neural_db1};	%for convenience
neural_bin_width = (h.skip + 1) / (h.speed_units / h.speed) * 1000;  %bin width of neural signal in ms

h2 = data.htb_header{neural_db2};	%for convenience
neural_bin_width2 = (h2.skip + 1) / (h2.speed_units / h2.speed) * 1000;  %bin width of neural signal in ms

AnalInterval = min( (StopEventBin + StopOffset) - (StartEventBin+ StartOffset) )  ;

low_freq_cutoff1 = 0;
high_freq_cutoff1 = 30; %hz low pass cut off for integrating power

low_freq_cutoff2 = 30;
high_freq_cutoff2 = 80; %hz low pass cut off for integrating power

low_freq_cutoff3 = 20;
high_freq_cutoff3 = 80; %hz low pass cut off for integrating power

%triangle correction factor for finite trial length
%triangle should be defined by length of spike train actually used for XCorr analysis
% triangle = (size(data.spike_data,2)*neural_bin_width/1000) - abs(time)/1000;
%spike_density_filter = [1/16 4/16 6/16 4/16 1/16];    
%spike_density_filter = [0.25 0.5 0.25];    
spike_density_filter = [1];
%trial = 0;

[conditions, unique_conds, param, num_conditions] = RegenerateConditionList(data, BegTrial, EndTrial, PATH, FILE, Protocol);


%now, remove trials from hor_disp and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:size(conditions,2);												% a vector of trial indices
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );  

num_rows = num_conditions + 2;
num_columns = 3;
plot_offset = num_columns * 2;

Fs = 1000;
NW = 5;
max_shift_predictor = 4;

output = 1;

curr_fig = figure;
set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [0 20 500 600], 'Name', 'Coherence Spectrogram');
subplot(10, 1, 1);

axis([0 100 0 100]);                    
axis('off');
xpos = -10;
ypos = 180;
font_size = 8;
bump_size = 30;

temp = strcat(PATH, FILE);
temp(temp == '\') = '/';
line = sprintf('Filename: %s %s', temp);
text(xpos,ypos, line,'FontSize',9);		ypos = ypos - bump_size;
line = sprintf('Spike Density Filter: %5.3f %5.3f %5.3f %5.3f %5.3f', spike_density_filter);
text(xpos,ypos, line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('Freq Band1 For Peak Power: %2d-%2d Hz', low_freq_cutoff1, high_freq_cutoff1 );
text(xpos,ypos, line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('Freq Band2 For Peak Power: %2d-%2d Hz', low_freq_cutoff2, high_freq_cutoff2 );
text(xpos,ypos, line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('Freq Band3 For Peak Power: %2d-%2d Hz', low_freq_cutoff3, high_freq_cutoff3 );
text(xpos,ypos, line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('Window Width %3d (ms)', segment_length );
text(xpos,ypos, line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('Window Step Size %3d (ms)', step_increment );
text(xpos,ypos, line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('Number of Tapers %d', NW * 2 - 1 );
text(xpos,ypos, line,'FontSize',font_size);		ypos = ypos - bump_size;


xpos = 40;
ypos = 150;
line = sprintf('Protocol: %d (%s)', Protocol, protocol_names{Protocol+1});
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('Trial range for analysis: %d -> %d', BegTrial, EndTrial);
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('Spike channels analyzed: %d, %d', SpikeChanA, SpikeChanB);
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('Analysis Start Event: %s, offset %d ms', event_names{start_code}, StartOffset);
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('Analysis Stop Event: %s, offset %d ms', event_names{stop_code}, StopOffset);
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;

for cond = 1: num_conditions
    if plot_PSTH == 1
        subplot(num_rows, num_columns, cond*num_columns + plot_offset - 2);
        phandle = plot(hist_time', hist_data2(cond,:), 'k-', hist_time', hist_data2B(cond,:), 'r-' );
        hold on;
        HIST_SCALE = 1.25 * max([max(max(hist_data2)) max(max(hist_data2B )) ]);
        xx = [mini_segment_start_times{cond}' mini_segment_start_times{cond}'];
        yy = [0 HIST_SCALE];
        if ~isempty(xx)
            plot(xx',yy', 'k--')
        end
        
        hold on;
        
        xx = [mini_segment_end_times{cond}' mini_segment_end_times{cond}'];
        yy = [0 max([max(hist_data2(cond,:)) max(hist_data2B(cond,:) ) ]) ];
        if ~isempty(xx)
            plot(xx',yy', 'g-.')
        end
        
        XLim([hist_time(1) hist_time(end)]);
        YLim([0 HIST_SCALE]);
    end
end

Selected_Average_Coherence = zeros(segment_length/2 + 1,num_conditions);
Coherency_Spectrogram = zeros(segment_length/2 + 1, length(segment_start_times), num_conditions );
Shuffle_Spectrogram = zeros(segment_length/2 + 1, length(segment_start_times), num_conditions );
Subtracted_Spectrogram = zeros(segment_length/2 + 1, length(segment_start_times), num_conditions );
Phase_Spectrogram = zeros(segment_length/2 + 1, length(segment_start_times), num_conditions );
Selected_Phase1 = zeros(1, num_conditions);
Selected_Phase2 = zeros(1, num_conditions);
Selected_Phase3 = zeros(1, num_conditions);
Selected_Coherence1 = zeros(1, num_conditions);
Selected_Coherence2 = zeros(1, num_conditions);
Selected_Coherence3 = zeros(1, num_conditions);

%analyze by unique condition
[E,V] = dpss(segment_length, NW);

for cond = 1 : num_conditions
    cond
    pars{cond} = NaN;  
    segment_ct = 0;
    Mean_Coherency_Vectors =  zeros(segment_length/2 + 1, 1);
    
    for segment = 1: length(segment_start_times)
        % These next few lines check for trials belonging to a single condition
        segment
        
        SetTrials(1,:) = (conditions(1,:) == unique_conds(cond,1) & select_trials); 	    
        for modality = 1: size(conditions,1) - 1         	
            NextSetTrials = (conditions(modality,:) == unique_conds(cond,modality) & select_trials); 
            SetTrials = SetTrials & NextSetTrials;				
        end	 
        
        reps = find(SetTrials==1); 
        num_reps = length(reps);		   
        
        start_spikebins = ceil( StartEventBin(reps) / neural_bin_width ) + ceil ( segment_start_times(segment)  / neural_bin_width ) - ceil (0.5 * segment_length/neural_bin_width);
        stop_spikebins = ceil( StartEventBin(reps) / neural_bin_width ) +  ceil ( segment_end_times(segment) / neural_bin_width  )- ceil (0.5 * segment_length/neural_bin_width);
        start_spikebins2 = ceil( StartEventBin(reps) / neural_bin_width2 ) + ceil ( segment_start_times(segment)  / neural_bin_width2 ) - ceil (0.5 * segment_length/neural_bin_width); 
        stop_spikebins2 = ceil( StartEventBin(reps) / neural_bin_width2 ) + ceil ( segment_end_times(segment) / neural_bin_width2  )- ceil (0.5 * segment_length/neural_bin_width);
        
        spikes = zeros(num_reps, segment_length);    
        spikes2 = zeros(num_reps, segment_length);
        Sxx_sum = 0;
        Syy_sum = 0;
        Sxy_sum = 0;

        % Now do cross correlograms 
        for trial = 1:num_reps           
            switch neural_db1
                case SPIKE_DB
                    spikes(trial,:) = conv(spike_density_filter, data.spike_data(SpikeChan, start_spikebins(trial): stop_spikebins(trial), reps(trial))  ); 
                    
                case LFP_DB
                    % use interp - spline function but constrained to go
                    % through data points
                    spikes(trial,:) = interp(data.lfp_data(SpikeChan,  start_spikebins(trial): stop_spikebins(trial), reps(trial)), neural_bin_width); 
                    spikes(trial,:) = detrend(spikes(trial,:));
            end   
            switch neural_db2
                case SPIKE_DB
                    spikes2(trial,:) = conv(spike_density_filter, data.spike_data(SpikeChan2, start_spikebins2(trial): stop_spikebins2(trial), reps(trial))  );     
                case LFP_DB
                    spikes2(trial,:) = interp(data.lfp_data(SpikeChan2, start_spikebins2(trial): stop_spikebins2(trial), reps(trial)), neural_bin_width2); 
                    spikes2(trial,:) = detrend(spikes2(trial,:));
            end   
        end
                
         mini_PSTH1 = mean(spikes,1);
         mini_PSTH2 = mean(spikes2,1);
        
        for trial = 1:num_reps           
              spikes(trial,:) = spikes(trial,:) - mini_PSTH1;
              spikes2(trial,:) = spikes2(trial,:) - mini_PSTH2;
            [Sxx, Syy, Sxy, freq] = multitaperCont(spikes(trial,:), spikes2(trial,:),E,V,[],Fs, 'unity', 'onesided');      
 %           [Sxx, Syy, Sxy, freq] = multitaper(spikes(trial,:), spikes2(trial,:),E,V,[],Fs, 'unity', 'onesided');      
            Sxx_sum = Sxx_sum + Sxx;
            Syy_sum = Syy_sum + Syy;
            Sxy_sum = Sxy_sum + Sxy;       
        end    
        
        Sxx_sum = Sxx_sum/num_reps;
        Syy_sum = Syy_sum/num_reps;
        Sxy_sum = Sxy_sum/num_reps;
        
        Coherency_Spectrogram(:, segment, cond) = (abs(Sxy_sum ./ ((Sxx_sum .* Syy_sum).^0.5 + 0.00000001) )).^2; 
        Phase_Spectrogram(:, segment, cond) = angle(Sxy_sum ./ ((Sxx_sum .* Syy_sum).^0.5 + 0.00000001) ); 
        
        Sxx2_sum = 0;
        Syy2_sum = 0;
        Sxy2_sum = 0;
        
        for shuffle = 1:max_shift_predictor
            shuffle_index = [num_reps - shuffle + 1: num_reps     1:num_reps - shuffle];
            for trial = 1:num_reps
                [Sxx, Syy, Sxy, freq] = multitaperCont(   spikes(trial,:), spikes2(shuffle_index(trial) ,:)  ,E,V,[], Fs, 'unity', 'onesided');      
%                [Sxx, Syy, Sxy, freq] = multitaper(   spikes(trial,:), spikes2(shuffle_index(trial) ,:)  ,E,V,[], Fs, 'unity', 'onesided');      
                Sxx2_sum = Sxx2_sum + Sxx;
                Syy2_sum = Syy2_sum + Syy;
                Sxy2_sum = Sxy2_sum + Sxy;       
            end           
        end
        
        Sxx2_sum = Sxx2_sum/(num_reps*max_shift_predictor );
        Syy2_sum = Syy2_sum/(num_reps*max_shift_predictor );
        Sxy2_sum = Sxy2_sum/(num_reps*max_shift_predictor );
        
        Shuffle_Spectrogram(:, segment, cond) = (abs(Sxy2_sum ./ ((Sxx2_sum .* Syy2_sum).^0.5 + 0.00000001) )).^2;   
        Subtracted_Spectrogram(:, segment,cond ) = (abs(   (Sxy_sum ./ ((Sxx_sum .* Syy_sum).^0.5 + 0.00000001)) - Shuffle_Spectrogram(:, segment, cond)   ) ).^2;    
        AutoSpectra1(:,cond) = Sxx_sum;
        AutoSpectra2(:,cond) = Syy_sum;
        for i = 1: length(mini_segment_start_times{cond})
            if (spec_time(segment) >= mini_segment_start_times{cond}(i) ) & (spec_time(segment) < mini_segment_end_times{cond}(i) )
%                Mean_Coherency_Vectors = Mean_Coherency_Vectors + Sxy_sum ./ ((Sxx_sum .* Syy_sum).^0.5 + 0.00000001) - Sxy2_sum ./ ((Sxx2_sum .* Syy2_sum).^0.5 + 0.00000001);    
                Mean_Coherency_Vectors = Mean_Coherency_Vectors + Sxy_sum ./ ((Sxx_sum .* Syy_sum).^0.5 + 0.00000001);    
                Selected_Average_Coherence(:,cond) = Selected_Average_Coherence(:,cond) + Subtracted_Spectrogram(:, segment, cond);
                segment_ct = segment_ct + 1;
            end
        end    
    end %end mini-segment
    
    Selected_Average_Coherence(:,cond) = Selected_Average_Coherence(:,cond) / segment_ct;
    
    Mean_Coherency_Vectors = Mean_Coherency_Vectors / segment_ct;
    Selected_Phase1(cond) = 180 / pi * angle( sum(Mean_Coherency_Vectors(logical(freq > low_freq_cutoff1) & logical(freq <= high_freq_cutoff1))) );
    Selected_Phase2(cond) = 180 / pi *angle( sum(Mean_Coherency_Vectors(logical(freq > low_freq_cutoff2) & logical(freq <= high_freq_cutoff2))) );
    Selected_Phase3(cond) = 180 / pi *angle( sum(Mean_Coherency_Vectors(logical(freq > low_freq_cutoff3) & logical(freq <= high_freq_cutoff3))) );
    Selected_Coherence1(cond) = abs( mean(Mean_Coherency_Vectors(logical(freq > low_freq_cutoff1) & logical(freq <= high_freq_cutoff1))) );
    Selected_Coherence2(cond) = abs( mean(Mean_Coherency_Vectors(logical(freq > low_freq_cutoff2) & logical(freq <= high_freq_cutoff2))) );
    Selected_Coherence3(cond) = abs( mean(Mean_Coherency_Vectors(logical(freq > low_freq_cutoff3) & logical(freq <= high_freq_cutoff3))) );
    
    Coherence = abs( Mean_Coherency_Vectors);
    Gamma_Peak_Freq(cond) = freq(find(Coherence == max(Coherence(freq >=20 & freq <=80 ))) );
    Gamma_Peak_Coherence(cond) = mean ( abs( Mean_Coherency_Vectors((freq >= Gamma_Peak_Freq(cond) - 5 & freq <= Gamma_Peak_Freq(cond) + 5) ) ) );
    
    num_selected_segments(cond) = length(mini_segment_start_times{cond});
end	% of analyzing for each condition

% plotting loop
for cond = 1: num_conditions
    subplot(num_rows, num_columns, cond*num_columns - 1 + plot_offset);
    phandle = plot(freq(freq <=200), Selected_Average_Coherence( (freq <= 200) ,cond) );
    xlabel ('Frequency (Hz)')
    if (cond == ceil(num_conditions /2)  )
        ylabel ('Coherence')
    end    
    if (cond == 1  )
        title ('Selected Segments')
    end    
    
    if (neural_db1 == LFP_DB ) && (neural_db2 == LFP_DB )
        xlim([0 100]);
        ylim([0 1]);
    else
        xlim([0 200]);
        ylim([0 max(max(Selected_Average_Coherence) )]);
    end   
    
    subplot(num_rows, num_columns, cond*num_columns + plot_offset);
    phandle = contourf(((1:size(Coherency_Spectrogram,2)) - 1)*step_increment, freq(freq <= 200), (Subtracted_Spectrogram( (freq <=200) ,:,cond)    ) );
    xlabel ('Time (ms)')
    if (cond == ceil(num_conditions /2)  )
        ylabel ('Frequency (Hz)')
    end    
    if (cond == 1  )
        title ('Shuffle subtracted')
    end    
    if (neural_db1 == LFP_DB ) && (neural_db2 == LFP_DB )
        ylim([0 100]);
    else
        ylim([0 200]);
    end
    colorbar  
    
end


%output tuning curve metrics
if (output == 1)
    i = size(PATH,2) - 1;
    while PATH(i) ~='\'	%Analysis directory is one branch below Raw Data Dir
        i = i - 1;
    end   
    PATHOUT = [PATH(1:i) 'Analysis\Coherence_Spectrograms\'];
    i = size(FILE,2) - 1;
    while FILE(i) ~='.'
        i = i - 1;
    end
    FILEOUT = [FILE(1:i - 1) '_' num2str(SpikeChanA) num2str(SpikeChanB) '.coh'];
    eval(['save ' PATHOUT FILEOUT  ' AutoSpectra1 AutoSpectra2 Gamma_Peak_Coherence Gamma_Peak_Freq max_shift_predictor Phase_Spectrogram num_selected_segments NW Subtracted_Spectrogram Selected_Coherence3 Selected_Coherence2 Selected_Coherence1 Selected_Phase1 Selected_Phase2 Selected_Phase3 Selected_Average_Coherence freq low_freq_cutoff1 high_freq_cutoff1 low_freq_cutoff2 high_freq_cutoff2 low_freq_cutoff3 high_freq_cutoff3 SpikeChanA SpikeChanB start_code stop_code BegTrial EndTrial StartOffset StopOffset StartEventBin StopEventBin PATH FILE'])
    
    %save figure
    FILEOUT = [FILE(1:i - 1) '_' num2str(SpikeChanA) num2str(SpikeChanB) '.fig'];        
    saveas (curr_fig, [PATHOUT FILEOUT]);
    
end