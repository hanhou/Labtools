%*******************************************************************************************************************
% BindingPermutationTest - Generates Cross Correlograms Smoothed Shuffled Correlograms, and subtracted 
%       correlograms.  Then does permutation test between hardcoded conditions.
% BJP 1/1/02

function BindingPermutationTest(data, Protocol, Analysis, SpikeChan, SpikeChan2, start_code, stop_code, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE);

line_types = ['b-'; 'r-'; 'g-'; 'k-'; 'g-'; 'g-'; 'g-'; 'g-'];
TEMPO_Defs;
ProtocolDefs;

%keep track of spike channels before reindexing
SpikeChanA = SpikeChan;
SpikeChanB = SpikeChan2;

switch SpikeChan
    case 4 %first LFP channel
        neural_db1 = LFP_DB;
        SpikeChan = 1;
    case 5 %second LFP channel
        neural_db1 = LFP_DB;
        SpikeChan = 2;
    case 1 %first Spike Channel
        neural_db1 = SPIKE_DB;
        SpikeChan = 1;
    case 2 %second Spike Channel
        neural_db1 = SPIKE_DB;
        SpikeChan = 2;
end

switch SpikeChan2
    case 4 %first LFP channel
        neural_db2 = LFP_DB;
        SpikeChan = 2;
    case 5 %second LFP channel
        neural_db2 = LFP_DB;
        SpikeChan = 1;
    case 1 %first Spike Channel
        neural_db2 = SPIKE_DB;
        SpikeChan = 2;
    case 2 %second Spike Channel
        neural_db2 = SPIKE_DB;
        SpikeChan = 1;
end
   
%first, we'll need the bin_width (in sec) of our spike raster + events log
h = data.htb_header{neural_db1};	%for convenience
neural_bin_width = (h.skip + 1) / (h.speed_units / h.speed) * 1000;  %bin width of neural signal in ms

AnalInterval = min( (StopEventBin + StopOffset) - (StartEventBin+ StartOffset) )  ;
  
rand('state',sum(100*clock));

range = 128;		% Range for generating cross correlograms
time = -range:+range; %actual time range for cross correlograms

peak_range = 2; % +/- range for width of averaging across peak
peak_window = 128;   % +/- range around 0 lag time in which peak is sought
high_freq_cutoff = 100;  %Hz lowpass cutoff for integrating power
low_freq_cutoff = 10;  %Hz highpass cutoff for integrating power

%triangle correction factor for finite trial length
triangle = (size(data.spike_data,2)*neural_bin_width/1000) - abs(time)/1000;
% spike_density_filter = [0.1 0.2 0.4 0.2 0.1];    
spike_density_filter = [0.25 0.5 0.25];    
% spike_density_filter = [1];

trial = 0;
  
num_modality = 0;
if SpikeChan == 1;
    SpikeChan2 = 2;   
else
    SpikeChan2 = 1;
end
    
[conditions, unique_conds, param, num_conditions] = RegenerateConditionList(data, BegTrial, EndTrial, PATH, FILE, Protocol);

cross_corr_peak_time = zeros(1,num_conditions);
cross_corr_peak_power = zeros(1,num_conditions);
sync_p_values = zeros(1,num_conditions);

%now, remove trials from hor_disp and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:size(conditions,2);												% a vector of trial indices
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );  
conditions = conditions(: , logical(select_trials));  
smoothed_shuffle = zeros(num_conditions, 2*range + 1);

overall_trial = 0;


if (Protocol == FIX_1_23)
    CondComp_1 = [2]; 
    CondComp_2 = [5];
end

if (Protocol == FIX_1_23_45)
    CondComp_1 = [2 2]; 
    CondComp_2 = [5 8];
end

if (Protocol == FIX_VARY_BACKCOLOR)
    CondComp_1 = [3 3 2]; 
    CondComp_2 = [2 5 5];
end

if (Protocol == FIX_VARY_HISTORY)
    CondComp_1 = [2 9]; 
    CondComp_2 = [5 10];
end
num_columns = 4;
num_rows = length(CondComp_1) + 1;
plot_offset = num_columns * 1;

num_bootstraps = 10000;
bootstrap_peak_power = zeros(num_bootstraps, num_conditions) ;

num_permutations = 10000;

output = 1;

figure   
set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [0 20 500 600], 'Name', 'Compare Correlation Levels');    
subplot(10, 1, 1);    
axis([0 100 0 100]);                    
axis('off');
xpos = 0;
ypos = 180;
font_size = 8;
bump_size = 30;

temp = strcat(PATH, FILE);
temp(temp == '\') = '/';
line = sprintf('Filename: %s %s', temp);
text(xpos,ypos, line,'FontSize',9);		ypos = ypos - bump_size;
% line = sprintf('Number of Repetitions: %3d', (size(data.spike_data,3)/num_conditions) );
% text(xpos,ypos, line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('Spike Density Filter: %5.3f %5.3f %5.3f %5.3f %5.3f', spike_density_filter);
text(xpos,ypos, line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('Window of Analysis: +/-%3d ms', peak_window);
text(xpos,ypos, line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('Range for Averaging Peak: %2d ms', 2*peak_range + 1);
text(xpos,ypos, line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('Freq Band For Integrating Power: %2d-%2d Hz', low_freq_cutoff, high_freq_cutoff );
text(xpos,ypos, line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('Number of Bootstraps: %8d', num_bootstraps);
text(xpos,ypos, line,'FontSize',font_size);		ypos = ypos - bump_size;

line = sprintf('Number of Permutations: %8d', num_permutations);
text(xpos,ypos, line,'FontSize',font_size);		ypos = ypos - bump_size;

%analyze by unique condition
for cond = 1 : num_conditions
    cond
    % These next few lines check for trials belonging to a single condition
    SetTrials(1,:) = (conditions(1,:) == unique_conds(cond,1) ); 	    
    for modality = 1: size(conditions,1) - 1         	
        NextSetTrials = (conditions(modality,:) == unique_conds(cond,modality) ); 
      	SetTrials = SetTrials & NextSetTrials;				
    end	 
     
    reps = find(SetTrials==1); 
    num_reps(cond) = length(reps);		   
    histo(cond, :) = zeros(1, size(time,2));
    histo2(cond,:) = zeros(1, size(time,2));
      
    spikes = [];    
    spikes2 = [];
 	% Now do cross correlograms 
    for trial = 1:num_reps(cond)
        start_spikebin = StartEventBin(reps(trial)) + StartOffset;
        switch neural_db1
            case SPIKE_DB
                spikes(trial,:) = conv(spike_density_filter, data.spike_data(SpikeChan, start_spikebin: start_spikebin + AnalInterval, reps(trial))  ); 
            case LFP_DB
                spikes(trial,:) = data.lfp_data(SpikeChan, floor(start_spikebin/neural_bin_width): floor( (start_spikebin)/neural_bin_width ) + floor( (AnalInterval)/neural_bin_width )  , reps(trial)); 
        end   
        switch neural_db2
            case SPIKE_DB
                spikes2(trial,:) = conv(spike_density_filter, data.spike_data(SpikeChan2, start_spikebin: start_spikebin + AnalInterval, reps(trial))  );     
            case LFP_DB
                spikes2(trial,:) = data.lfp_data(SpikeChan2, floor(start_spikebin/neural_bin_width): floor( (start_spikebin)/neural_bin_width ) + floor( (AnalInterval)/neural_bin_width )   , reps(trial)); 
        end   
    end   
  
    corr_matrix = zeros( num_reps(cond)*(num_reps(cond)), 2*range + 1);

    
    combin = 0;
    for trial1 = 1 : num_reps(cond)
        for trial2 = 1 : num_reps(cond)
            combin = combin + 1;
            corr_matrix(combin,:) = xcorr(spikes(trial1, :), spikes2(trial2, : ), range)./triangle;		    
            %normalize by geometric mean spike rate
            norm_factor = ( (mean( spikes(trial1, :) ) *1000  ) * (mean( spikes2(trial2, : ) )*1000 ) )^0.5;
            if (norm_factor ~= 0)
                corr_matrix(combin,:) = corr_matrix(combin,:) / norm_factor;
            end    

        end
    end       
    
    % first obtain the general structure of the correlogram by creating smoothed shuffle
    % only from trial pairings that aren't cross correlograms
    crosses = [];
    crosses(linspace(1,(num_reps(cond) )^2, num_reps(cond) )) = 1;
    smoothed_shuffle(cond,:) = mean(corr_matrix(~logical(crosses) , :)  );
        
    for trial1 = 1: size(corr_matrix,1)
        %first subtract shuffle
        corr_matrix(trial1,:) = corr_matrix(trial1,:) - smoothed_shuffle(cond,:);
        %subtract DC
        corr_matrix(trial1,:) = corr_matrix(trial1,:) - mean(corr_matrix(trial1,:) );
    end
        
    cross_corr_matrix{cond}(:,:) = corr_matrix(logical(crosses),:);
    histo(cond,:) = mean(cross_corr_matrix{cond}(:,:) ); 

    cross_corr_peak_time(cond) = max(-1 + find(time == -peak_window) + time(find( histo(cond,find(time == -peak_window):find(time == peak_window) ) == max(histo(cond, find(time == -peak_window):find(time == peak_window) ) ) )  ) );
    [freq2, ampl2] = FourierTransform_1D(0.001*[1 : (size(histo,2) - 1)], histo(cond, 2:end), size(histo,2) - 1, 1, 0);
    cross_corr_peak_power(cond) = sum(ampl2( (freq2 > low_freq_cutoff) & (freq2 < high_freq_cutoff) )  .* ampl2( (freq2 > low_freq_cutoff) & (freq2 < high_freq_cutoff) )  );
  
    ampl5 = zeros(1, (size(corr_matrix,2) - 1)/2 );    
      
    for strap = 1:num_bootstraps
        shuffle = zeros(1, size(time,2));
        %shuffle trials   with replacement
        shuffle_index = fix(rand(1,num_reps(cond) )*combin) + 1;
        shuffle = mean(corr_matrix(shuffle_index,:)  );
        bootstrap_peak_time(strap, cond) = max( time(-1 + find(time == -peak_window) + find( shuffle(find(time == -peak_window):find(time == peak_window)) == max(shuffle(find(time == -peak_window):find(time == peak_window) ) ) ) )  );
        [freq, ampl] = FourierTransform_1D(0.001*[1 : (length(shuffle) - 1)], shuffle(2:end), length(shuffle) - 1, 1, 0);   
        bootstrap_peak_power(strap, cond) = sum(ampl((freq < high_freq_cutoff) & (freq > low_freq_cutoff)  ) .* ampl( (freq < high_freq_cutoff)& (freq > low_freq_cutoff)   )   );

        ampl5 = ampl5 + ampl.^2;       
    end  

    ampl5 = ampl5 / num_bootstraps;
    smooth_corr_peak_power(cond) = sum(ampl5(  (freq2 < high_freq_cutoff) & (freq > low_freq_cutoff) )   );
    smooth_corr_peak_power(cond) = 0;
    
    cross_corr_peak_power(cond) = cross_corr_peak_power(cond) - smooth_corr_peak_power(cond);
    bootstrap_peak_power(:,cond) = bootstrap_peak_power(:,cond) - smooth_corr_peak_power(cond);    
    
    sync_p_values(cond) = mean(bootstrap_peak_power(:, cond) > cross_corr_peak_power(cond) );
  
  
end	% of analyzing for each condition
 
comp_p_values = zeros(1,length(CondComp_1) );
perm_peak_diff = zeros(num_permutations, length(CondComp_1) );
    
%should also preallocate for permcorrpeak power
perm_corr_peak_power = zeros(num_permutations, 2);
perm_corr_peak_time = zeros(num_permutations, 2);

corr_max = 1.5*max(max(histo) );
corr_min = 1.5*min(min(histo) );

for comp = 1: length(CondComp_1) 
    corr_peak_diff(comp) = cross_corr_peak_power(CondComp_1(comp) ) - cross_corr_peak_power(CondComp_2(comp) );
    trials_cond1 = ones(1, num_reps(CondComp_1(comp)) );
    trials_cond2 = ones(1, num_reps(CondComp_2(comp)) );
    for perm = 1 : num_permutations
        comp_matrix1 = 0*cross_corr_matrix{CondComp_1(comp)};
        comp_matrix2 = 0*cross_corr_matrix{CondComp_2(comp)};         
        trials_cond1(1:min([num_reps(CondComp_1(comp))  num_reps(CondComp_2(comp)) ]) ) =  fix(2 * rand(1, min([num_reps(CondComp_1(comp))  num_reps(CondComp_2(comp)) ]) ) );
        trials_cond2(1:min([num_reps(CondComp_1(comp))  num_reps(CondComp_2(comp)) ]) )= trials_cond1(1:min([num_reps(CondComp_1(comp))  num_reps(CondComp_2(comp)) ]) );

        comp_matrix1(1:sum(trials_cond1) , :) = cross_corr_matrix{CondComp_1(comp)}(logical(trials_cond1),:);
        comp_matrix1(sum(trials_cond1) + 1 : size(comp_matrix1,1) , :) = cross_corr_matrix{CondComp_2(comp)}(logical(~trials_cond2),:);
        comp_matrix2(1:sum(trials_cond2) , :) = cross_corr_matrix{CondComp_2(comp)}(logical(trials_cond2),:);
        comp_matrix2(sum(trials_cond2) + 1 : size(comp_matrix2,1) , :) = cross_corr_matrix{CondComp_1(comp)}(logical(~trials_cond1),:);
        
        %subtract DC offset
        correlogram1 = mean(comp_matrix1) - mean(mean(comp_matrix1) );
        correlogram2 = mean(comp_matrix2) - mean(mean(comp_matrix2) );
        perm_corr_peak_time(perm, 1) = max(-1 + find(time == -peak_window) + time(find( correlogram1(find(time == -peak_window):find(time == peak_window) ) == max(correlogram1(find(time == -peak_window):find(time == peak_window) ) ) )  ) );
%        perm_corr_peak_heights(perm, 1) = mean(correlogram1(find(time == perm_corr_peak_time(perm,1) ) - peak_range : find(time == perm_corr_peak_time(perm,1) ) + peak_range   ),2)';
        %calculate area
%        peak = histo(cond,find(time == perm_corr_peak_time(perm, 1) ) - peak_range : find(time == perm_corr_peak_time(perm, 1) ) + peak_range   );
%        perm_corr_peak_heights(perm, 1) = sum(peak(peak > 0 ) );
 
        %calculate power in correlogram
        [freq, ampl] = FourierTransform_1D(0.001*[1 : (length(correlogram1) - 1)], correlogram1(2:end), length(correlogram1) - 1, 1, 0);   
        perm_corr_peak_power(perm, 1) = sum(ampl((freq < high_freq_cutoff) & (freq > low_freq_cutoff)  ) .*ampl( (freq < high_freq_cutoff) & (freq > low_freq_cutoff)  ) );

        perm_corr_peak_time(perm, 2) = max(-1 + find(time == -peak_window) + time(find( correlogram2(find(time == -peak_window):find(time == peak_window) ) == max(correlogram2(find(time == -peak_window):find(time == peak_window) ) ) )  ) );
%        perm_corr_peak_heights(perm, 2) = mean(correlogram2(find(time == perm_corr_peak_time(perm,2) ) - peak_range : find(time == perm_corr_peak_time(perm,2) ) + peak_range   ),2)';
%        peak = histo(cond,find(time == perm_corr_peak_time(perm, 2) ) - peak_range : find(time == perm_corr_peak_time(perm, 2) ) + peak_range   );
%        perm_corr_peak_heights(perm, 2) = sum(peak(peak > 0 ) );

        %calculate power in correlogram
        [freq, ampl3] = FourierTransform_1D(0.001*[1 : (length(correlogram2) - 1)], correlogram2(2:end), length(correlogram2) - 1, 1, 0);   
        perm_corr_peak_power(perm, 2) = sum(ampl3((freq < high_freq_cutoff) & (freq > low_freq_cutoff) ) .*ampl3(  (freq < high_freq_cutoff) & (freq > low_freq_cutoff) )  );

        perm_peak_diff(perm, comp) = perm_corr_peak_power(perm, 1) - perm_corr_peak_power(perm, 2);
        if (perm < 100)
            subplot(num_rows, num_columns, 4*comp - 3 + plot_offset)
            plot(time, correlogram1, 'g' );
            hold on;
            
            subplot(num_rows, num_columns,4*comp - 2 + plot_offset)
            plot(time, correlogram2, 'g' );           
            hold on;
        end
    end    
    comp_p_values(comp) = mean( perm_peak_diff(:,comp) > corr_peak_diff(comp) );
%     mean_diffs(comp) = mean(perm_peak_diff(:,comp) );
%     histo2 = zeros(num_permutations,2);
%     histo2(:,1) = perm_peak_diff(:,comp);
%     histo2(:,2) = NaN;
%     histo2(1:num_permutations*0.05, 2) = corr_peak_diff(comp);
    
    subplot(num_rows, num_columns, 4*comp - 3 + plot_offset)
    plot(time,histo(CondComp_1(comp), :), 'b' );
    xlim([min(time) max(time)] );
    ylim([corr_min corr_max] );
    
    subplot(num_rows, num_columns, 4*comp - 2 + plot_offset)
    plot(time,histo(CondComp_2(comp), :), 'k');
    xlim([min(time) max(time)] );
    ylim([corr_min corr_max] );
 
    subplot(num_rows, num_columns, 4*comp - 1 + plot_offset);
    axis([0 100 0 100]);                    
    axis('off');
    xpos = -20;
    ypos = 100;
    font_size = 8;
    bump_size = 25;
    line = '';
    for modality = 1: length(param) 
        line = [line param{modality}, num2str(unique_conds( CondComp_1(comp),modality)), ' '];
    end   
    text(xpos,ypos, line,'FontSize',font_size);		ypos = ypos - bump_size;

    line = sprintf('Bootstrap P Value Cond 1: %7.5f', sync_p_values( CondComp_1(comp) ));
    text(xpos,ypos, line, 'FontSize',font_size);		ypos = ypos - bump_size;
  
    line = '';
    for modality = 1: length(param) 
        line = [line param{modality}, num2str(unique_conds( CondComp_2(comp),modality)), ' '];
    end   
    text(xpos,ypos, line,'FontSize',font_size);		ypos = ypos - bump_size;
    
    line = sprintf('Bootstrap P Value Cond 1: %7.5f', sync_p_values( CondComp_2(comp) ));
    text(xpos,ypos, line, 'FontSize',font_size);		ypos = ypos - bump_size;
    
    line = sprintf('Permutation P Value: %7.5f', comp_p_values(comp));
    text(xpos,ypos, line, 'FontSize',font_size);		ypos = ypos - bump_size;
    axis('off');
    hold on;

end


    %output tuning curve metrics
    if (output == 1)
        i = size(PATH,2) - 1;
        while PATH(i) ~='\'	%Analysis directory is one branch below Raw Data Dir
            i = i - 1;
        end   
        PATHOUT = [PATH(1:i) 'Analysis\Correlograms\'];
        i = size(FILE,2) - 1;
        while FILE(i) ~='.'
            i = i - 1;
        end
        FILEOUT = [FILE(1:i) 'prm'];
        
        eval(['save ' PATHOUT FILEOUT  ' cross_corr_peak_power cross_corr_peak_time bootstrap_peak_power bootstrap_peak_time sync_p_values comp_p_values num_bootstraps perm_corr_peak_power perm_corr_peak_time num_permutations corr_peak_diff perm_peak_diff range peak_range peak_window spike_density_filter low_freq_cutoff high_freq_cutoff SpikeChanA SpikeChanB neural_db1 neural_db2 start_code stop_code BegTrial EndTrial StartOffset StopOffset StartEventBin StopEventBin PATH FILE'])
    end

[corr_peak_diff; comp_p_values]

