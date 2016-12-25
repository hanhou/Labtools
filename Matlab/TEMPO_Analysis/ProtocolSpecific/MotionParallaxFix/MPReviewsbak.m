%-----------------------------------------------------------------------------------------------------------------------
%-- MPReviews.m -- Analyses to placate reviewers on the Nature paper.  Looks at the extra signal added to the RM 
% response to create the MP response (MP-RM) for each depth including null.  Returns amplitude and phase and 
% amplitude significance using MPBootstrap3.  Also returns the difference in phase for each depth relative to the
% null.  In addition, monotonicity results for all conditions are tacked into this script.       
%-- Started by JWN, 12/05/07
%-- Last by JWN, 12/05/07
%-----------------------------------------------------------------------------------------------------------------------
function MPReviewsbak(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE);

ver = '1.1';
TEMPO_Defs;
Path_Defs;
symbols = {'bo' 'rs' 'gd' 'kv' 'm<' 'c>' 'bv' 'rv'};
line_types2 = {'b--' 'r--' 'g--' 'k--' 'g.-' 'b.-' 'r-.' 'k.'};
line_types4 = {'b-' 'r-' 'g-' 'k-' 'm-' 'c-' 'y-' 'b-'};
line_types5 = {'bo-' 'rs-' 'gd-' 'kv-' 'm<-' 'c>-' 'yo-' 'bs-'};
NULL_VALUE = -9999;

disp(sprintf('(MPReviews v%s) Started at %s.',ver,datestr(now,14)));

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
    disp('(MPReviews) Fatal Error: Two phases required to calculate modulation indices.');
    return;
end
trials = size(MPphase,2);

% Get the mean firing rates for all the trials
area = 'MST';  % Kluge! 80 for MT and 80 for MST (see Kruse et al 2002), +80 for transfer function delay
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
modulation_phase = zeros(num_trial_types, num_depths)+888;
modulation_sig = zeros(num_trial_types, num_depths)+888;
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
%         modulation_fft = fft(saved_counts(i,j,:));
%         modulation_amp(i,j) = abs(modulation_fft(2)); 
%         modulation_phase(i,j) = 180+angle(modulation_fft(2))*180/pi;
%         modulation_sig(i,j) = MPBootstrap3(saved_counts(i,j,:), modulation_amp(i,j));
    end
end
% Now begin subtractions and FFT calculations
% First MP-RM
MPRMcounts = saved_counts(1,:,:)-saved_counts(3,:,:);
for j=1:num_depths  % Ten MPdepths (which includes null)
    MPRM_fft = fft(MPRMcounts(1,j,:));
    MPRM_amp(j) = abs(MPRM_fft(2));
    MPRM_phase(j) = 180+angle(MPRM_fft(2))*180/pi;
    MPRM_sigp(j) = (MPBootstrap3(MPRMcounts(1,j,:), MPRM_amp(j)))<0.025;
end
% Normalize MPRM_amps from whatever they are in to be relative to the max MPRM_amp
MPRM_amp = MPRM_amp./max(MPRM_amp);
% Get phase differences for a histogram
for j=1:num_depths-1
    phases = unwrap([MPRM_phase(1)/(180/pi) MPRM_phase(j+1)/(180/pi)]);
    phase_difference(j) = abs(phases(1)-phases(2))*180/pi;  %in degrees
end

% Test TCs for monotonicity (3 ways TC can be judged monotonic)
% numparamq = 3; % Polynomial fit, used in mono3Reviews
numparamq = 2; % Binomial fit, to see if results are different with that
numparaml = 1; % Linear fit (1 or 2 fewer parameters)
range = [-2:.5:2];
monotonic = zeros(1,6)+888;
for i = 1:6
    if(isempty(find(uMPtrial_types==i-1))) continue;  end;  % Break out if none from that condition
    reps = floor(sum(MPtrial_types == i-1)/(num_depths*num_phase));  % Moved reps in here because different conditions may now have different numbers of reps.
    mean_data = zeros(reps*2,num_depths-1);
    for j = 1:num_depths-1
        tmp = spike_rates(MPdepths == uMPdepths(j+1) & MPtrial_types == i-1)';
        mean_data(:,j) = tmp(1:reps*2);  % ignore extra incomplete reps
    end
    a = mean(mean_data);
    if(issorted(a)|issorted(a(9:-1:1)))
        monotonic(i) = 2;  % TC is strictly monotonic returns "2"
    else
        pf = polyfit(range,a,numparamq);
        q = polyval(pf,range);
        if(issorted(q)|issorted(q(9:-1:1)))
            monotonic(i) = 1;  % 2. Polynomial fit is strictly monotonic returns "1"
        else
            pf1 = polyfit(range,a,numparaml);
            l = polyval(pf1,range);
            sseq = sum((a - q).^2);
            ssel = sum((a - l).^2);
            F = ((ssel-sseq)/(numparamq-numparaml))/(sseq/(9-numparamq));
            monotonic(i) = 1 - fcdf(F,(numparamq-numparaml),(9-numparamq));  % 3. P>0.05 means polynomial fit not better than linear, returns P
        end
    end
end

% Write results for this cell to 1 file
PATHOUT = 'Z:\Data\MOOG\Ovid\Analysis\';
filenames = {'Reviews'};
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
        fprintf(fid, 'paMPRMnull paMPRMn20 paMPRMn15 paMPRMn10 paMPRMn05 paMPRM00 paMPRM05 paMPRM10 paMPRM15 paMPRM20 ');
        fprintf(fid, 'sMPRMnull sMPRMn20 sMPRMn15 sMPRMn10 sMPRMn05 sMPRM00 sMPRM05 sMPRM10 sMPRM15 sMPRM20 ');
        fprintf(fid, 'sumsigps ');
        fprintf(fid, 'pdn20 pdn15 pdn10 pdn05 pd00 pd05 pd10 pd15 pd20 ');
        fprintf(fid, 'monoMP monoBD monoRM monoC monoEO monoHO');        
        fprintf(fid, '\r\n');
    end
    fprintf(fid,'%10s', strtok(FILE,'.'));
    fprintf(fid,' %+2.5f', monkid, cellid, MPRM_phase, MPRM_sigp, sum(MPRM_sigp), phase_difference, monotonic);
    fprintf(fid,'\r\n');
    fprintf(fid,'%10s', strtok(FILE,'.'));
    fprintf(fid,' %+2.5f', monkid, cellid, MPRM_amp, MPRM_sigp);
    fprintf(fid,'\r\n');
    fclose(fid);
end
% final_data = [bins squeeze(MPRMcounts)'];
% PATHOUT = 'Z:\Data\MOOG\Ovid\Analysis\';
% filenames = {'MPRMPSTHs'};
% for i = 1:1
%     outfile = cell2mat(strcat(PATHOUT,strtok(FILE,'.'),'_',filenames(i),'.txt'));
%     fid = fopen(outfile, 'w');  % Open text file.
%     fprintf(fid, 'Bins Null neg2.0 neg1.5 neg1.0 neg0.5 0.0 0.5 1.0 1.5 2.0');
%     fprintf(fid, '\r\n');
%     for j = 1:size(final_data,1)
%         fprintf(fid,' %+2.4f', final_data(j,:));
%         fprintf(fid,'\r\n');
%     end
%     fclose(fid);
% end
disp('(MPReviews) Done.');
return;