%-----------------------------------------------------------------------------------------------------------------------
%-- MPDirSelect.m -- Analyses to placate reviewers on the Nature paper
%-- Started by JWN, 12/21/07
%-- Last by JWN, 12/21/07
%-----------------------------------------------------------------------------------------------------------------------
function MPDirSelect(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE);

ver = '1.0';
TEMPO_Defs;
Path_Defs;
symbols = {'bo' 'rs' 'gd' 'kv' 'm<' 'c>' 'bv' 'rv'};
line_types2 = {'b--' 'r--' 'g--' 'k--' 'g.-' 'b.-' 'r-.' 'k.'};
line_types4 = {'b-' 'r-' 'g-' 'k-' 'm-' 'c-' 'y-' 'b-'};
line_types5 = {'bo-' 'rs-' 'gd-' 'kv-' 'm<-' 'c>-' 'yo-' 'bs-'};
NULL_VALUE = -9999;

disp(sprintf('(MPDirSelect v%s) Started at %s.',ver,datestr(now,14)));

[monkid, cellid, runstr]=strread(FILE,'m%dc%dr%s.htb');
% Get the trial type, depth values, and movement phase for each condition in the condition_list[]
MPdepths = data.moog_params(PATCH_DEPTH,:,MOOG);
uMPdepths = unique(MPdepths);
num_depths = size(uMPdepths,2);
MPtrial_types = data.moog_params(MP_TRIAL_TYPE,:,MOOG);
uMPtrial_types = unique(MPtrial_types);  % Conditions present
%Place breakouts here.  This is what SelectiveAnalysis could be all about!
%if(isempty(find(uMPtrial_types==###))) return;  end;
%if(isempty(find(uMPtrial_types==0))) disp('(MPSelectiveAnalysis) Breakout: No MP');  return;  end;  % BREAKOUT ENABLED!

num_trial_types = length(uMPtrial_types);
MPphase = data.moog_params(MOVEMENT_PHASE,:,MOOG);
uMPphase = unique(MPphase);
num_phase = size(uMPphase,2);
if(num_phase ~= 2)
    disp('(MPDirSelect) Fatal Error: Two phases required to calculate modulation indices.');
    return;
end
trials = size(MPphase,2);

% Get the mean firing rates for all the trials
area = 'MT';  % Kluge! 80 for MT and 80 for MST (see Kruse et al 2002), +80 for transfer function delay
if(strcmp(area,'MT'))  % Don't change this one!
    latency = 160;  % MT guess
else
    latency = 160;  % MST guess
end 
begin_time = find(data.event_data(1,:,1)==StartCode) + latency; % Each trial always has the same start time so may as well use trial 1
end_time = begin_time + 1999; % 2s trial
if(max(max(max(data.spike_data))) > 1)
    data.spike_data = cast(data.spike_data>0,'double');
end
raw_spikes = data.spike_data(1,begin_time:end_time,:);
spont_spikes = data.spike_data(1,begin_time-500:begin_time,:);
spike_rates = 1000*squeeze(mean(raw_spikes))';  % The hard way
spont_rates = 1000*squeeze(mean(spont_spikes))';
total_spike_bins = end_time - begin_time;
num_reduced_bins = 39;
bin_width = total_spike_bins/(num_reduced_bins+1);  % ~2000ms/(39+1) = ~50ms;

% Look at modulations
% Get differenced PSTHs for each condition and save them (saved_counts), along with modulation amplitudes, phases, and significance
saved_counts = zeros(num_trial_types, num_depths, num_reduced_bins+1);
modulation_amp = zeros(num_trial_types, num_depths)+888;
modulation_DC = zeros(num_trial_types, num_depths)+888;
% modulation_sig = zeros(num_trial_types, num_depths)+888;
for i = [1 3]
    if(isempty(find(uMPtrial_types==i-1))) continue;  end;  % Break out if none from that condition
    for j=1:num_depths  % Ten MPdepths (which includes null)
        indices0 = logical((MPtrial_types == i-1) & (MPdepths == uMPdepths(j)) & (MPphase == uMPphase(1)));
        indices180 = logical((MPtrial_types == i-1) & (MPdepths == uMPdepths(j)) & (MPphase == uMPphase(2)));
        raw_spikes = data.spike_data(1,begin_time:end_time,indices0);
        hist_data = sum(raw_spikes,3);
        [bins, counts0] = SpikeBinner(hist_data, 1, bin_width, 0);
        counts0 = counts0*(1000/bin_width)/sum(indices0);  % Convert to instantaneous firing rates
        raw_spikes = data.spike_data(1,begin_time:end_time,indices180);
        hist_data = sum(raw_spikes,3);
        [bins, counts180] = SpikeBinner(hist_data, 1, bin_width, 0);
        counts180 = counts180*(1000/bin_width)/sum(indices180);  % Convert to instantaneous firing rates
        % Now subtract phases
        saved_counts(i,j,:) = counts0-counts180;
        modulation_fft = fft(saved_counts(i,j,:));
        modulation_amp(i,j) = abs(modulation_fft(2))/length(bins)/0.5; % half amplitude in spikes/s
        modulation_DC(i,j) = mean([counts0' counts180']); % in spikes/s
        % modulation_sig(i,j) = MPBootstrap3(saved_counts(i,j,:), abs(modulation_fft(2)));
    end
end
F1DC = modulation_amp./modulation_DC;
F1DCMP = F1DC(1,:);
F1DCRM = F1DC(3,:);
final_data = reshape([F1DCMP;F1DCRM],[1,20]);

% Write results for this cell to 1 file
PATHOUT = 'Z:\Data\MOOG\Ovid\Analysis\';
filenames = {'DirSelect'};
for i = 1:1
    outfile = cell2mat(strcat(PATHOUT,area,'_',filenames(i),'.txt'));
    headerflag = 0;
    if (exist(outfile) == 0) % File does not yet exist, so print a header
        headerflag = 1;
    end
    fid = fopen(outfile, 'a');  % Open text file.
    if (headerflag)
        fprintf(fid, 'FILE ');
        fprintf(fid, 'monkid cellid ');
        fprintf(fid, 'F1DCMPnull F1DCRMnull F1DCMPn20 F1DCRMn20 F1DCMPn15 F1DCRMn15 F1DCMPn10 F1DCRMn10 F1DCMPn05 F1DCRMn05 F1DCMP00 F1DCRM00 F1DCMP05 F1DCRM05 F1DCMP10 F1DCRM10 F1DCMP15 F1DCRM15 F1DCMP20 F1DCRM20 ');
        fprintf(fid, 'F1DCMPavg F1DCRMavg ');  % Includes null!
        fprintf(fid, 'absMPiPDI absRMiPDI');  % Not calculated here; need to c&p these in later  
        fprintf(fid, '\r\n');
    end
    fprintf(fid,'%10s', strtok(FILE,'.'));
    fprintf(fid,' %+2.5f', monkid, cellid, final_data, mean(F1DCMP), mean(F1DCRM));
    fprintf(fid,'\r\n');
    fclose(fid);
end
disp('(MPDirSelect) Done.');
return;