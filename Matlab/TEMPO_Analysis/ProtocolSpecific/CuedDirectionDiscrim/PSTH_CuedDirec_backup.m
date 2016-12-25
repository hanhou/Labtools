%-----------------------------------------------------------------------------------------------------------------------
%-- PSTH_CuedDirec.m -- Plot PSTHs for each stimulus condition
%--	VR, 9/21/05
%-----------------------------------------------------------------------------------------------------------------------

function PSTH_CuedDirec(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

TEMPO_Defs;		
Path_Defs;
ProtocolDefs;	%needed for all protocol specific functions - contains keywords - BJP 1/4/01

%get the column of values of directions in the dots_params matrix
direction = data.dots_params(DOTS_DIREC,:,PATCH1);
unique_direction = munique(direction');
Pref_direction = data.one_time_params(PREFERRED_DIRECTION);

%get the motion coherences
coherence = data.dots_params(DOTS_COHER, :, PATCH1);
unique_coherence = munique(coherence');

%get the cue validity: -1=Invalid; 0=Neutral; 1=Valid; 2=CueOnly
cue_val = data.cue_params(CUE_VALIDITY,:,PATCH2);
unique_cue_val = munique(cue_val');
cue_val_names = {'NoCue','Invalid','Neutral','Valid','CueOnly'};
NOCUE = -2; INVALID = -1; NEUTRAL = 0; VALID = 1; CUEONLY = 2;

%get the cue directions
cue_direc = data.cue_params(CUE_DIREC, :, PATCH1);
unique_cue_direc = munique(cue_direc');
%cue_dir_type = 1 if PrefDir, 0 if Neutral Cue, -1 if Null Cue
cue_dir_type = logical( (squeeze_angle(cue_direc) == Pref_direction) & (cue_val ~= NEUTRAL) ) - logical( (squeeze_angle(cue_direc) ~= Pref_direction) & (cue_val ~= NEUTRAL) );
unique_cue_dir_type = munique(cue_dir_type');
cue_dir_typenames = {'Null','Neutral','Pref'};

%compute cue types - 0=neutral, 1=directional, 2=cue_only
cue_type = abs(cue_val); %note that both invalid(-1) and valid(+1) are directional
unique_cue_type = munique(cue_type');

%get indices of any NULL conditions (for measuring spontaneous activity)
null_trials = logical( (coherence == data.one_time_params(NULL_VALUE)) );

% keyboard

%now, select trials that fall between BegTrial and EndTrial
trials = 1:length(coherence);
%a vector of trial indices
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );

%get outcome for each trial: 0=incorrect, 1=correct
trial_outcomes = logical (data.misc_params(OUTCOME,:) == CORRECT);
trial_choices = ~xor((direction==Pref_direction),trial_outcomes); %0 for Null Choices, 1 for Pref Choices

linetypes = {'b-','r-','g-'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%first align the psths with the stimulus events - cue onset, motion onset

%first compute the psth centered around the cue onset, and a psth centered around the stimulus onset
precue = 200; %time to show before the cue starts
postcue = 400; %time to show after cue starts
prestim = 300; %time to display before the visual stimulus
poststim = 1100; %time to display after the visual stimulus
binwidth = 25; %in ms (used for psth)
bw = 50; %in ms, used for roc
spksamprate = 1000;
stddev = 15; %in ms, std dev of guassian used for filtering
buff = 3*stddev; 
gaussfilt = normpdf([1:2*buff+1],buff+1,stddev); %gaussian filter 3 std.dev's wide
timing_offset = 71; %in ms, time between VSTIM_ON_CD and detectable light on screen; use to offset motion onset.


for k = 1:length(unique_cue_dir_type)
    for i = 1:length(unique_coherence)
        for j = 1:2 %PrefChoice=1, NullChoice=2
            %first select the relevant trials and get a raster
            select = trials(select_trials & (coherence == unique_coherence(i)) & ...
                (trial_choices == 2-j) & (cue_dir_type == unique_cue_dir_type(k)) );
            for m = 1:length(select)
                t_stimon = find(data.event_data(1,:,select(m)) == VSTIM_ON_CD) + timing_offset;
                prestim_raster{i,j,k}(m,:) = data.spike_data(SpikeChan, t_stimon-prestim-buff:t_stimon+poststim+buff, select(m));
            end
            prestim_psth{i,j,k} = sum(prestim_raster{i,j,k},1)./length(select).*spksamprate; %psth is NOT binned
            sm_prestim_psth{i,j,k} = conv(gaussfilt, prestim_psth{i,j,k}); %convolve with the gaussian filter
            sm_prestim_psth{i,j,k} = sm_prestim_psth{i,j,k}(2*buff+1:end-2*buff); %lop off the edges
            
%             % THIS IS OLD CODE FOR BINNING
%             prestimraster_length = length(prestim_raster{i,j,k,1});
%             trunc_prestimraster_length = prestimraster_length - mod(prestimraster_length,binwidth);
%             prestim_psth{i,j,k}=zeros(1,trunc_prestimraster_length/binwidth);
%             %below is where the rasters are binned, summed and averaged.
%             for m = 1:length(select)
%                 binned_prestim_raster{i,j,k}(m,:) = sum(reshape(prestim_raster{i,j,k,m}(1:trunc_prestimraster_length),binwidth,trunc_prestimraster_length/binwidth));
%                 prestim_psth{i,j,k} = prestim_psth{i,j,k} + binned_prestim_raster{i,j,k}(m,:);
%             end
%             prestim_psth{i,j,k} = prestim_psth{i,j,k}./length(select)./binwidth.*spksamprate;
        end
        %now compute a running roc metric for the two choices
        for v = 1:size(sm_prestim_psth{i,1,k},2) %time bin
            if isempty(prestim_raster{i,1,k}) | isempty(prestim_raster{i,2,k})
                prestim_roc{i,k}(v) = NaN;
            else
                t_start = v + buff+1 - bw/2-1; %first time point of bin of width bw centered around time v
                t_stop = v + buff+1 + bw/2;
                pc = sum(prestim_raster{i,1,k}(:,t_start:t_stop),2); %spikes at the bin centered around vth time point for pref choice trials
                nc = sum(prestim_raster{i,2,k}(:,t_start:t_stop),2); %spikes at the bin centered around vth time point for null choice trials
                prestim_roc{i,k}(v) = rocn(pc,nc,100);
            end
        end
    end
    %repeat this (w/o roc) collapsing across all coherences for the cue response
    for j = 1:2 %again for the two choices - a little kludgey organization but allows roc computation
        select = trials(select_trials & (trial_choices == 2-j) & (cue_dir_type == unique_cue_dir_type(k)) );
        for m = 1:length(select)
            t_cueon = find(data.event_data(1,:,select(m)) == CUE_ON_CD) + timing_offset;
            postcue_raster{j,k}(m,:) = data.spike_data(SpikeChan, t_cueon-precue-buff:t_cueon+postcue+buff, select(m));
        end
        postcue_psth{j,k} = sum(postcue_raster{j,k},1)./length(select).*spksamprate; %psth is NOT binned
        sm_postcue_psth{j,k} = conv(gaussfilt, postcue_psth{j,k}); %convolve with the gaussian filter
        sm_postcue_psth{j,k} = sm_postcue_psth{j,k}(2*buff+1:end-2*buff); %lop off the edges
        
%         %THIS IS OLD STUFF FOR BINNING
%         postcueraster_length = length(postcue_raster{j,k,1});
%         trunc_postcueraster_length = postcueraster_length - mod(postcueraster_length,binwidth);
%         postcue_psth{j,k}=zeros(1,trunc_postcueraster_length/binwidth);
%         for m = 1:length(select)
%             binned_postcue_raster{j,k}(m,:) = sum(reshape(postcue_raster{j,k,m}(1:trunc_postcueraster_length),binwidth,trunc_postcueraster_length/binwidth));
%             postcue_psth{j,k} = postcue_psth{j,k} + binned_postcue_raster{j,k}(m,:);
%         end
%         postcue_psth{j,k} = postcue_psth{j,k}./length(select)./binwidth.*spksamprate;
    end    
end

%find max of the psths
maxy_prestim = zeros(length(unique_coherence),1);
for j = 1:length(unique_coherence)
    temp = sm_prestim_psth(j,:,:);
    for i = 1:prod(size(temp))
        if (max(temp{i})>maxy_prestim(j))
            maxy_prestim(j) = max(temp{i});
        end
    end
end

%now plot the peristimulus psths
postcue_x = [-precue:(length(sm_postcue_psth{1,1,1})-precue-1)];
prestim_x = [-prestim:(length(sm_prestim_psth{1,1,1})-prestim-1)];
h(1)=figure;
set(h(1),'PaperPosition', [.2 .2 8 10.7], 'Position', [250 50 500 573], 'Name', sprintf('%s: Peri-Stimulus Time Histogram',FILE));

for j = 1:2 %PrefChoice = 1; NullChioce=2;
    %first plot the peri-CueOnset psth, then underneath plot the peri-VStim psths
    subplot(1+length(unique_coherence), 2, j); hold on;
    for k = 1:length(unique_cue_dir_type)
        plot(postcue_x,sm_postcue_psth{j,k},linetypes{k});
    end
    axis tight;
    xlabel('Time about Cue Onset');
    if j==1
        ylabel('F.R.{Hz}');
        title(sprintf('%s: PrefDir Choices',FILE));
%         legh=legend('NullDir','Neutral','PrefDir','Location','NorthEast');
%         set(legh,'box','off');
    else
        title('NullDir Choices');
    end
    for i = 1:length(unique_coherence)
        subplot(1+length(unique_coherence), 2, i*2+j);
        hold on;
        for k = 1:length(unique_cue_dir_type)
            if ~( isempty(sm_prestim_psth{i,1,k}) | isempty(sm_prestim_psth{i,2,k}) )
                plot(prestim_x,sm_prestim_psth{i,j,k},linetypes{k});
            end
        end
        axis tight;
        ylim([0 maxy_prestim(i)]);
        if j==1
            ylabel(sprintf('Coh= %3.1f%%',unique_coherence(i)));
        end
        if i==1
            
        elseif i==length(unique_coherence)
            xlabel('Time about VStim Onset (ms)');
        end
    end
end

%now plot the roc time courses per coherence
h(2)=figure;
set(h(2),'PaperPosition', [.2 .2 8 10.7], 'Position', [250 50 500 573], 'Name', sprintf('%s: ROC',FILE));
for i = 1:length(unique_coherence)
    subplot(length(unique_coherence),1,i); hold on;    
    for k = 1:length(unique_cue_dir_type)
        plot(prestim_x,prestim_roc{i,k},linetypes{k});
    end
    axis tight
    plot(xlim, [0.5 0.5], 'k:');
    if i==1
        title(sprintf('%s: ROC values, sorted by cue direction',FILE));
    elseif i==length(unique_coherence)
        xlabel('Time about VStim Onset (ms)');
    end
    ylabel(sprintf('Coh = %6.1f',unique_coherence(i)));
    ylim([0 1]);
end

% keyboard
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%now repeat this around saccades
%first compute the psth centered around the saccade onset
presacc = 400; %time to show before the saccade starts
postsacc = 200; %time to show after saccade starts


for k = 1:length(unique_cue_dir_type)
    for i = 1:length(unique_coherence)
        for j = 1:2 %PrefChoice=1, NullChoice=2
            %first select the relevant trials and get a raster
            select = trials(select_trials & (coherence == unique_coherence(i)) & ...
                (trial_choices == 2-j) & (cue_dir_type == unique_cue_dir_type(k)) );
            for m = 1:length(select)
                t_sacc = find(data.event_data(1,:,select(m)) == SACCADE_BEGIN_CD);
                sacc_raster{i,j,k}(m,:) = data.spike_data(SpikeChan, t_sacc-presacc-buff:t_sacc+postsacc+buff, select(m));
            end
            
            sacc_psth{i,j,k} = sum(sacc_raster{i,j,k},1)./length(select).*spksamprate; %psth is NOT binned
            sm_sacc_psth{i,j,k} = conv(gaussfilt, sacc_psth{i,j,k}); %convolve with the gaussian filter
            sm_sacc_psth{i,j,k} = sm_sacc_psth{i,j,k}(2*buff+1:end-2*buff); %lop off the edges

%             %OLD STUFF FOR BINNING
%             saccraster_length = length(sacc_raster{i,j,k,1});
%             trunc_saccraster_length = saccraster_length - mod(saccraster_length,binwidth);
%             sacc_psth{i,j,k}=zeros(1,trunc_saccraster_length/binwidth);
%             %below is where the rasters are binned, summed and averaged.
%             for m = 1:length(select)
%                 binned_sacc_raster{i,j,k}(m,:) = sum(reshape(sacc_raster{i,j,k,m}(1:trunc_saccraster_length),binwidth,trunc_saccraster_length/binwidth));
%                 sacc_psth{i,j,k} = sacc_psth{i,j,k} + binned_sacc_raster{i,j,k}(m,:);
%             end
%             sacc_psth{i,j,k} = sacc_psth{i,j,k}./length(select)./binwidth.*spksamprate;
        end

        %now compute a running roc metric for the two choices
        for v = 1:size(sm_sacc_psth{i,1,k},2)
            t_start = v + buff+1 - bw/2-1; %first time point of bin of width bw centered around time v
            t_stop = v + buff+1 + bw/2;
            pc = sum(sacc_raster{i,1,k}(:,t_start:t_stop),2); %spikes at the bin centered around vth time point for pref choice trials
            nc = sum(sacc_raster{i,2,k}(:,t_start:t_stop),2); %spikes at the bin centered around vth time point for null choice trials
            sacc_roc{i,k}(v) = rocn(pc,nc,100);
        end
    end
end
%find max of the psths
maxy = zeros(length(unique_coherence),1);
for j = 1:length(unique_coherence)
    temp = sm_sacc_psth(j,:,:);
    for i = 1:prod(size(temp))
        if (max(temp{i})>maxy(j))
            maxy(j) = max(temp{i});
        end
    end
end
%now plot the peristimulus psths
sacc_x = [-presacc:(length(sm_sacc_psth{1,1})-presacc-1)];
h(3)=figure;
set(h(3),'PaperPosition', [.2 .2 8 10.7], 'Position', [250 50 500 573], 'Name', sprintf('%s: Peri-Saccadic Time Histogram',FILE));
for j = 1:2 %PrefChoice = 1; NullChioce=2;
    for i = 1:length(unique_coherence)
        subplot(length(unique_coherence), 2, (i-1)*2+j);
        hold on;
        for k = 1:length(unique_cue_dir_type)
            plot(sacc_x,sm_sacc_psth{i,j,k},linetypes{k});
        end
        axis tight;
        ylim([0 maxy(i)])
        plot([0 0],ylim,'k');
        if j==1
            ylabel(sprintf('Coh= %3.1f%%',unique_coherence(i)));
        end
        if i==1
            if j==1
                title(sprintf('%s: PrefDir Choices',FILE));
            else
                title('NullDir Choices');
            end
        elseif i==length(unique_coherence)
            xlabel('Time about Saccade Onset (ms)');
        end
    end
end

%now plot the roc time courses per coherence
h(2)=figure;
set(h(2),'PaperPosition', [.2 .2 8 10.7], 'Position', [250 50 500 573], 'Name', sprintf('%s: Saccade-aligned ROC',FILE));
for i = 1:length(unique_coherence)
    subplot(length(unique_coherence),1,i); hold on;    
    for k = 1:length(unique_cue_dir_type)
        plot(sacc_x,sacc_roc{i,k},linetypes{k});
    end
    axis tight
    plot(xlim, [0.5 0.5], 'k:');
    if i==1
        title(sprintf('%s: ROC values, sorted by cue direction',FILE));
    elseif i==length(unique_coherence)
        xlabel('Time about Saccade Onset (ms)');
    end
    ylabel(sprintf('Coh = %6.1f',unique_coherence(i)));
    ylim([0 1]);
end

end



