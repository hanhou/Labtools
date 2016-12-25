function SimulateSpikeTrains(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset,  StartEventBin, StopEventBin, PATH, FILE);   

%%%%%%%%%%%%%%%%%%
%%%%%%%%%%  Neuron Parameter
%%%%%%%%%%%%%%%%%%%
anal_bin_width = 1; %ms
corr_jitter =2; % ms lambda for gaussian jitter for sync pulses
corr_lag = 0; % ms time lag between coordinated spikes
%spike_density_filter = [0.25 0.5 0.25];
spike_density_filter = [1];
oscil_freq = 50;

gen_PSTH = 2;
warning off MATLAB:divideByZero;

if gen_PSTH == 1
    PSTH_bin_width = 40;
    %% PSTH contains function normalized to 1 for spike rate
    sin_freq = [1 1]; %Hz
    dc_offset = [0.5 0.5];
    amplitude = [1 1];
    phase = [-pi/2 -pi/2];
    trial_length = 3; %seconds
    num_reps = 30;
    FR = [100 100];
    
    data_hist_time = PSTH_bin_width/2 + 0:PSTH_bin_width: trial_length*1000;
    %data segments centered on maxima of PSTH

    data_mini_segment_start_times = [350 450 550 1350 1450 1550 2350 2450 2550];
    data_mini_segment_end_times = data_mini_segment_start_times + 100;
    PATHOUT = ['z:\data\Simulations\Model\' num2str(corr_jitter) '\'];

else
    %load PSTH
    i = size(PATH,2) - 1;
    while PATH(i) ~='\'	%Analysis directory is one branch below Raw Data Dir
        i = i - 1;
    end   
    PATHOUT = [PATH(1:i) 'Analysis\PSTH\'];
    i = size(FILE,2) - 1;
    while FILE(i) ~='.'
        i = i - 1;
    end
    FILEOUT = [FILE(1:i - 1)  '_12.psh'];
    eval (['load ' PATHOUT FILEOUT ' -MAT']); 
    
    %Load ISI Histograms
    i = size(PATH,2) - 1;
    while PATH(i) ~='\'	%Analysis directory is one branch below Raw Data Dir
        i = i - 1;
    end   
    PATHOUT = [PATH(1:i) 'Analysis\ISI\'];
    i = size(FILE,2) - 1;
    while FILE(i) ~='.'
        i = i - 1;
    end
    FILEOUT = [FILE(1:i - 1)  '.isi'];
    eval (['load ' PATHOUT FILEOUT ' -MAT']); 

    %reset filename
    i = size(PATH,2) - 1;
    while PATH(i) ~='\'	%Analysis directory is one branch below Raw Data Dir
        i = i - 1;
    end   
    while PATH(i) ~='\'	%Analysis directory is one branch below Raw Data Dir
        i = i - 1;
    end   
    PATHOUT = [PATH(1:i) 'Raw\'];
    i = size(FILE,2) - 1;
    while FILE(i) ~='.'
        i = i - 1;
    end
    FILEOUT = [FILE(1:i - 1)  '.htb'];
    i = size(PATHOUT,2) - 1;
    while PATHOUT(i) ~='\'	%Analysis directory is one branch below Raw Data Dir
        i = i - 1;
    end   
    PATHOUT = [PATHOUT(1:i-1) '\Analysis\Simulated_Spike_Data\'];
    i = size(FILEOUT,2) - 1;
    while FILEOUT(i) ~='.'
        i = i - 1;
    end
    FILEOUT = [FILEOUT(1:i - 1) '_s000'];
            FR = [0 0];
end

%% correlation parameters
peak_range = 5; % +/- range for width of averaging across peak
peak_window = 10;   % +/- range around 0 lag time in which peak is sought
corr_range = 64; %ms +/- bounds for correlograms
corr_parms = [peak_range peak_window corr_range];
high_freq_cutoff = 100;
low_freq_cutoff = 10;
num_iterations = 2;
time = -corr_range:anal_bin_width : corr_range;
gabor_fit = 0;

%bootstrap parameters
num_bootstraps = 10000;
num_iterations = 1;
corr_parms = [peak_range peak_window corr_range gabor_fit];
corr_levels = [0.0 0.03 0.05 0.1];

neuron_parms(2) = corr_jitter;
neuron_parms(3) = corr_lag;
neuron_parms(4) = num_iterations;
neuron_parms(5) = FR(1);
neuron_parms(6) = FR(2);
neuron_parms(7) = oscil_freq;

simulation = 1;
for sim = 1:length(corr_levels)

    neuron_parms(1) = corr_levels(sim);
    
    while exist([PATHOUT FILEOUT '.spk']);
        FILENUM = str2num(FILEOUT(end-2:end) );
        FILENUM = FILENUM + 1;
        FILE = num2str(FILENUM);
        FILEOUT(end-2:end)  = '000';
        FILEOUT(end - length(FILE) + 1:end) = FILE;
    end    
    
    if (gen_PSTH ~= 1)
        cond = 2;
        PSTH_bin_width = 40;
        data_ISI(1,1:301) = histc( isi{cond},0:300);
        data_ISI(2,1:301) = histc( isi2{cond},0:300);
        data_raster(1,:) = sum( hist_data{cond} );
        data_raster(2,:) = sum( hist_dataB{cond} );
        data_total_spikes(1,:) = total_spikes1(cond);
        data_total_spikes(2,:) = total_spikes2(cond);
        data_PSTH(1,:) = hist_data2(cond,:);
        data_PSTH(2,:) = hist_data2B(cond,:);
        data_mini_segment_end_times = mini_segment_end_times{cond};
        data_mini_segment_start_times = mini_segment_start_times{cond};
        data_hist_time = PSTH_bin_width/2 + 0:PSTH_bin_width: size(data_raster,2);

    end
     
    [PATHOUT FILEOUT]
    clear sparse_raster spike_raster bootstrap_peak_heights bootstrap_peak_times corr_matrix;
    %    [spike_raster, spike_stats, actual_corr_rate, sparse_raster]= gen_spike_trains(PSTH, neuron_parms, corr_rates, num_reps, trial_length, anal_bin_width, spike_density_filter, PATHOUT, FILEOUT);
    
    [spike_raster spike_raster2 sim_PSTH] = gen_spike_trains(neuron_parms, anal_bin_width, data_PSTH, data_ISI, data_total_spikes, data_raster, PATHOUT, FILEOUT);
    
    [NCC(sim)] = gen_corr_matrix(spike_raster, spike_raster2, anal_bin_width, data_hist_time, data_mini_segment_start_times, data_mini_segment_end_times, corr_parms, neuron_parms, PATHOUT, FILEOUT);
     [Coh1(sim), Coh2(sim), Coh3(sim)] = gen_spike_coherence(spike_raster, spike_raster2, sim_PSTH, data_hist_time, data_mini_segment_start_times, data_mini_segment_end_times, corr_parms, neuron_parms, PATHOUT, FILEOUT); 
    
     while exist([PATHOUT FILEOUT '.spk']);
        FILENUM = str2num(FILEOUT(end-2:end) );
        FILENUM = FILENUM + 1;
        FILE = num2str(FILENUM);
        FILEOUT(end-2:end)  = '000';
        FILEOUT(end - length(FILE) + 1:end) = FILE;
    end    
    
    if (gen_PSTH ~= 1)
        cond = 5;
        PSTH_bin_width = 40;
        data_ISI(1,1:301) = histc( isi{cond},0:300);
        data_ISI(2,1:301) = histc( isi2{cond},0:300);
        data_raster(1,:) = sum( hist_data{cond} );
        data_raster(2,:) = sum( hist_dataB{cond} );
        data_total_spikes(1,:) = total_spikes1(cond);
        data_total_spikes(2,:) = total_spikes2(cond);
        data_PSTH(1,:) = hist_data2(cond,:);
        data_PSTH(2,:) = hist_data2B(cond,:);
        data_mini_segment_end_times = mini_segment_end_times{cond};
        data_mini_segment_start_times = mini_segment_start_times{cond};
        data_hist_time = PSTH_bin_width/2 + 0:PSTH_bin_width: size(data_raster,2);

    end
     
    [PATHOUT FILEOUT]
    clear sparse_raster spike_raster bootstrap_peak_heights bootstrap_peak_times corr_matrix;
    %    [spike_raster, spike_stats, actual_corr_rate, sparse_raster]= gen_spike_trains(PSTH, neuron_parms, corr_rates, num_reps, trial_length, anal_bin_width, spike_density_filter, PATHOUT, FILEOUT);
    
    [spike_raster spike_raster2 sim_PSTH] = gen_spike_trains(neuron_parms, anal_bin_width, data_PSTH, data_ISI, data_total_spikes, data_raster, PATHOUT, FILEOUT);
    
    [NCC(sim)] = gen_corr_matrix(spike_raster, spike_raster2, anal_bin_width, data_hist_time, data_mini_segment_start_times, data_mini_segment_end_times, corr_parms, neuron_parms, PATHOUT, FILEOUT);
    [Coh1(sim), Coh2(sim), Coh3(sim)] = gen_spike_coherence(spike_raster, spike_raster2, sim_PSTH, data_hist_time, data_mini_segment_start_times, data_mini_segment_end_times, corr_parms, neuron_parms, PATHOUT, FILEOUT); 
       
    simulation = simulation + 1;
end    %sim


