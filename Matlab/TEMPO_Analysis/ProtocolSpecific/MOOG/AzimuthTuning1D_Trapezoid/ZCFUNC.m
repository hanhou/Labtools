function     ZCFUNC(data, Protocol, Analysis, SpikeChan, StartEventBin, StopEventBin,StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
TEMPO_Defs
Path_Defs;
ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP
DataMatrixDef;
Stand_D = 0.1;
trials = 1:size(data.moog_params,2);
select_trials = false(size(trials));
colorset = {[0.1 0.1 1 0.9],[1 0.1 0.1 0.9],[0 0.5 0 0.9]};
dt = 1/500;
if length(BegTrial) == 1 && BegTrial > 0 % Backward compatibility
    select_trials(BegTrial:EndTrial) = true;
elseif all(BegTrial > 0) % To be included
    select_trials(BegTrial) = true;
elseif all(BegTrial < 0) % To be excluded
    select_trials(-BegTrial) = true;
    select_trials = ~ select_trials;
else
    disp('Trial selection error...');
    keyboard;
end
% Time information
SPIKE_DB = 2; %ZC 2021125
h = data.htb_header{SPIKE_DB};
%h.skip + 1 sample num in one bin
spike_timeWin = 1000 * (h.skip + 1) / (h.speed_units / h.speed); % in ms % HH20140520

% Trial information
[temp_a,temp_t1,temp_t2,temp_v,temp_d,Mon_Outcome,RightAns,Mon_Choice,total_trials,One_target_trail,star_Denstity,star_Coherence,star_LifeTime,stim_type,sa,sas,Tcof] = ZC_Trail_Info(data);

distance_per_trial = temp_d(select_trials);
A_per_trial = temp_a(select_trials);
T1_per_trial = temp_t1(select_trials);
T2_per_trial = temp_t2(select_trials);
for i = 1:sum(select_trials)
    T1BinTime1 = floor(T1_per_trial(i)*500*(1-Tcof(i))); %500hz ZC 
    T1BinTime2 = floor(T1_per_trial(i)*500*(1+Tcof(i))); 
    T2BinTime = floor(T2_per_trial(i)*500); 
    A_profile_per_trail{i} = A_per_trial(i)*([ones(1,T1BinTime1)/(1-Tcof(i)),zeros(1,T2BinTime),-ones(1,T1BinTime2)/(1+Tcof(i))]);
    V_profile_per_trail{i} = cumsum(A_profile_per_trail{i})*dt;
    D_profile_per_trail{i} = cumsum(V_profile_per_trail{i})*dt;
end
Outcome_per_trail = Mon_Outcome(select_trials);
unique_distance = munique(distance_per_trial');
unique_distance(isnan(unique_distance))=[];
repetitionN = floor(length(distance_per_trial) /length(unique_distance)) ;


% Spike data
spike_in_bin = squeeze(data.spike_data(SpikeChan,:,select_trials));   % 5000 * TrialNum for default
spike_in_bin(spike_in_bin > 2) = 2; % Some MU from Chan21 would have more than one spike per timebin. HH20150210
% Event data
events_in_bin = squeeze(data.event_data(1,:,select_trials));
hold off
for i= 1:3
selected_trail = stim_type ==i;
[Aspike{i},Vspike{i},Dspike{i}] = get_par_alined_spike(A_profile_per_trail,V_profile_per_trail,D_profile_per_trail,events_in_bin,spike_in_bin,selected_trail);
%Align neron data with a,v,d across time 
figure(4)
plot(Aspike{i}(1,:),Aspike{i}(2,:),'Color',colorset{i})
hold on
title('align with A')
legend('VES','VIS','COM')
figure(5)
plot(Vspike{i}(1,:),Vspike{i}(2,:),'Color',colorset{i})
hold on
title('align with V')
legend('VES','VIS','COM')
figure(6)
plot(Dspike{i}(1,:),Dspike{i}(2,:),'Color',colorset{i})
hold on
title('align with D')
legend('VES','VIS','COM')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
align_markers = [ VSTIM_ON_CD, VSTIM_OFF_CD, SACCADE_BEGIN_CD];  % Desired markers: target onset & target offset & sac onset
psth_windows = [-500 2000; -1000 1500; -2000 500]; % Pre-trigger and Post-trigger PSTH windows (cover all possible range)
MarkerTitle = {'Stim On','Stim OFF','Sac On'};
Plot_Time_Win = {-500/2:2000/2,-1000/2:1500/2,-2000/2:500/2};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% Anne 2014 NN: 10 ms bin + 50 ms Gaussian smoothing. HH20141207
binSize = 10; % in ms
stepSize = 10; % in ms
smoothFactor = 50; % in ms (Gaussian kernel)

for j = 1:length(align_markers)    % For each desired marker
    [align_offsets(:,j),~] = find(events_in_bin == align_markers(j));  % Fast way to find offsets of each trial
    
    % --- New Smoothed PSTH for each trial. HH20141207 ---
    % Align spikes    
    spike_aligned{j} = nan(sum(select_trials),ceil((psth_windows(j,2) - psth_windows(j,1)) / spike_timeWin) + 1);
    
    for i = 1:sum(select_trials)
        winBeg = align_offsets(i,j) + ceil(psth_windows(j,1) / spike_timeWin);
        winEnd = align_offsets(i,j) + ceil(psth_windows(j,2) / spike_timeWin);
        if winBeg > 0 && winEnd <= size(spike_in_bin,1) % OK time windows
            spike_aligned{j}(i,:) = spike_in_bin (winBeg : winEnd,i)';   % Align each trial
        else
            beep;
            fprintf('WARNING: PSTH window error at spike_aligned{%g}(%g,:)\n',j,i);
%             error('PSTH window error at spike_aligned: Line 207');
%             winBeg
%             spike_aligned{j}(i, 2-winBeg : end) = spike_in_bin (1 : winEnd,i)';
%             spike_aligned{j}(i,:) = NaN;
        end
    end
end
    
%Long VS Short 
%Only correct trials.
LCidx = distance_per_trial  > Stand_D & Outcome_per_trail == 1;
SCidx = distance_per_trial <= Stand_D & Outcome_per_trail == 1;

for j = 1:length(align_markers)    % For each desired marker
    for modality = 1:3            
        modidx = stim_type == modality;
        SpikeLC{j,modality} = smooth(mean(spike_aligned{j}(LCidx & modidx,:)),50);
        SpikeSC{j,modality} = smooth(mean(spike_aligned{j}(SCidx & modidx,:)),50);
    end
end

%    spike_aligned{j}(:,:) = 
 %   MeanSpk_aligned{j}    
for j = 1:length(align_markers) 
    figure(j)
    hold off
    for modality = 1:3                
        plot(Plot_Time_Win{j},SpikeLC{j,modality},'-','LineWidth',1.5,'Color',colorset{modality})        
        hold on    
        plot(Plot_Time_Win{j},SpikeSC{j,modality},'--','LineWidth',1.5,'Color',colorset{modality})
        title(MarkerTitle{j})
    end
end


