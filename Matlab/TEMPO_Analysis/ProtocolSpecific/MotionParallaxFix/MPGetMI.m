%-----------------------------------------------------------------------------------------------------------------------
%-- MPGetMI.m -- Returns modulation indices for every phase pair. 
%-- Started by JWN, 3/2/05
%-- Last by JWN, 8/2/07
%-----------------------------------------------------------------------------------------------------------------------
function [mis] = MPGetMI(latency, data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, FILE);

TEMPO_Defs;		%needed for defines like IN_T1_WIN_CD
ProtocolDefs;	%needed for all protocol specific functions - contains keywords - BJP 1/4/01
Path_Defs;

disp(sprintf('(MPGetMI) Started at %s.',datestr(now,14)));
% Repeat much of MP_PSTH

% Get the trial type, depth values, and movement phase for each condition in the condition_list[]
MPdepths = data.moog_params(PATCH_DEPTH,:,MOOG);
MPtrial_types = data.moog_params(MP_TRIAL_TYPE,:,MOOG);
MPphase = data.moog_params(MOVEMENT_PHASE,:,MOOG);

begin_time = find(data.event_data(1,:,1)==StartCode) + latency;
% end_time = find(data.event_data(1,:,1)==StopCode) + latency;
end_time = begin_time + 1999; % 2s trial
if(max(max(max(data.spike_data))) > 1)
    data.spike_data = cast(data.spike_data>0,'double');
end
% raw_spikes = data.spike_data(1,begin_time:end_time,:);  % Each trial always has the same start time so may as well use trial 1
% spike_rates = 1000*squeeze(mean(raw_spikes))';  % The hard way
total_spike_bins = end_time - begin_time;
num_reduced_bins = 39;
bin_width = total_spike_bins/(num_reduced_bins+1);  % ~2000ms/(39+1) = ~50ms;

uMPtrial_types = unique(MPtrial_types);
uMPdepths = unique(MPdepths);
uMPphase = unique(MPphase);
num_trial_types = size(uMPtrial_types,2);
num_depths = size(uMPdepths,2);
num_phase = size(uMPphase,2);
if(num_phase ~= 2)
    disp('(MPGetMI) Fatal Error: Two phases required to calculate modulation indices.');
    return;
end
count_data = zeros(num_trial_types,num_depths,num_phase,num_reduced_bins+1);
subtracted_count_data = zeros(num_trial_types,num_depths,num_reduced_bins+1);
trials = size(MPphase,2);
reps = floor(trials/(num_trial_types*num_depths*num_phase));
mis = zeros(num_trial_types,num_depths,reps);
phases = zeros(num_trial_types,num_depths,reps);

% Loop through to get binned spikes, differences between phases, and then MIs
for i=1:num_trial_types  % Four MPtrial_type blocks (0=MP,1=BD,2=RM,3=C)
    for j=1:num_depths  % Ten MPdepths (which includes null)
        indices0 = logical((MPtrial_types == uMPtrial_types(i)) & (MPdepths == uMPdepths(j)) & (MPphase == uMPphase(1)));
        trials0 = find(indices0==1);
        indices180 = logical((MPtrial_types == uMPtrial_types(i)) & (MPdepths == uMPdepths(j)) & (MPphase == uMPphase(2)));
        trials180 = find(indices180==1);
        for k=1:reps    % Calculate MIs from paired trials
            % Get counts for phase 0
            raw_spikes = data.spike_data(1,begin_time:end_time,trials0(k));
            [bins, counts0] = SpikeBinner(raw_spikes, 1, bin_width, 0);
            % Get counts for phase 180
            raw_spikes = data.spike_data(1,begin_time:end_time,trials180(k));
            [bins, counts180] = SpikeBinner(raw_spikes, 1, bin_width, 0);
            % Subtract and calculate MI
            counts = counts180-counts0;
            phase_fft = fft(counts);
            phases(i,j,k) = phase_fft(2);
            full_fft = abs(phase_fft);  % Counts per half-cycle... 
            % mis(i,j,k) = full_fft(2);  % ...which for 2s = counts/second, so we could just do this, but we make Matlab do extra work...
            mis(i,j,k) = full_fft(2)/(size(full_fft,1)/2);  % Convert to counts/bin
            mis(i,j,k) = mis(i,j,k) * 1000/bin_width;  % Convert counts/bin to counts/second
        end
    end
end

% Only look at MP condition
% near_phase = angle(sum(sum((phases(1,2:5,:)))))*180/pi;
% near_conf = 1000/abs(sum(sum(phases(1,2:5,:))));
% far_phase = angle(sum(sum((phases(1,7:10,:)))))*180/pi;
% far_conf = 1000/abs(sum(sum(phases(1,7:10,:))));

% Write results for this cell to a common file
% area = 'MT';
% PATHOUT = 'Z:\Data\MOOG\Barracuda\Analysis\';
% filenames = {'PhaseData'};
% for i = 1:1
%     outfile = cell2mat(strcat(PATHOUT,area,'_',filenames(i),'.txt'));
%     headerflag = 0;
%     if (exist(outfile) == 0) % File does not yet exist, so print a header
%         headerflag = 1;
%     end
%     fid = fopen(outfile, 'a');  % Open text file.
%     if (headerflag)
%         fprintf(fid, 'FILE NearPhase NearConf FarPhase FarConf');
%         fprintf(fid, '\r\n');
%     end
%     fprintf(fid,'%10s', strtok(FILE,'.'));
%     fprintf(fid,' %+2.4f', near_phase, near_conf, far_phase, far_conf);
%     fprintf(fid,'\r\n');
%     fclose(fid);
% end
disp('(MPGetMI) Done.');
return;