%-----------------------------------------------------------------------------------------------------------------------
%-- MPFisher.m -- Calculates Fisher Information-based thresholds around 0.0 simulated depth for MP, BD, RM, and C conditions
% Could easily be expanded to add analysis of EO and HO.
%-- Started by JWN, 11/12/08
%-- Last by JWN, 11/12/08
%-----------------------------------------------------------------------------------------------------------------------
function MPFisher(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE);

ver = '1.0';
TEMPO_Defs;
Path_Defs;
symbols = {'bo' 'rs' 'gd' 'kv' 'm<' 'c>' 'bv' 'rv'};
line_types2 = {'b--' 'r--' 'g--' 'k--' 'g.-' 'b.-' 'r-.' 'k.'};
line_types4 = {'b-' 'r-' 'g-' 'k-' 'm-' 'c-' 'y-' 'b-'};
line_types5 = {'bo-' 'rs-' 'gd-' 'kv-' 'm<-' 'c>-' 'yo-' 'bs-'};
NULL_VALUE = -9999;

disp(sprintf('(MPFisher v%s) Started at %s.',ver,datestr(now,14)));

% Get the trial type, depth values, and movement phase for each condition in the condition_list[]
MPdepths = data.moog_params(PATCH_DEPTH,:,MOOG);
uMPdepths = unique(MPdepths);
num_depths = size(uMPdepths,2);
MPtrial_types = data.moog_params(MP_TRIAL_TYPE,:,MOOG);
uMPtrial_types = unique(MPtrial_types);  % Conditions present
%Place breakouts here.  This is what SelectiveAnalysis could be all about!
%if(isempty(find(uMPtrial_types==###))) return;  end;
if(isempty(find(uMPtrial_types==3))) disp('(MPFisher) Breakout: No C');  return;  end;  % BREAKOUT ENABLED!

num_trial_types = length(uMPtrial_types);
MPphase = data.moog_params(MOVEMENT_PHASE,:,MOOG);
uMPphase = unique(MPphase);
num_phase = size(uMPphase,2);
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
spike_rates = 1000*squeeze(mean(raw_spikes))';  % The hard way

for i = 1:4 % MP=1,BD=2,RM=3,C=4
   indices = logical((MPtrial_types == (i-1)) & (MPdepths == -0.5));
   near05 = mean(spike_rates(indices));
   indices = logical((MPtrial_types == (i-1)) & (MPdepths == 0.5));
   far05 = mean(spike_rates(indices));
   slopeat0(i) = near05 - far05;
   indices = logical((MPtrial_types == (i-1)) & (MPdepths == 0));
   varat0 = var(spike_rates(indices));
   FI(i) = slopeat0(i)^2/varat0;
   thresholdat0(i) = 1/sqrt(FI(i));
end

% Write results for this cell to 1 file
PATHOUT = 'Z:\Data\MOOG\Ovid\Analysis\';
filenames = {'FI'};
for i = 1:1
    outfile = cell2mat(strcat(PATHOUT,area,'_',filenames(i),'.txt'));
    headerflag = 0;
    if (exist(outfile) == 0) % File does not yet exist, so print a header
        headerflag = 1;
    end
    fid = fopen(outfile, 'a');  % Open text file.
    if (headerflag)
        fprintf(fid, 'FILE ');
        fprintf(fid, 'FIMP FIBD FIRM FIC threshMP threshBD threshRM threshC slopeMP slopeBD slopeRM slopeC ratCoMP ratCoBD congruency');
        fprintf(fid, '\r\n');
    end
    fprintf(fid,'%10s', strtok(FILE,'.'));
    fprintf(fid,' %+2.5f', FI, thresholdat0, slopeat0, thresholdat0(4)/thresholdat0(1), thresholdat0(4)/thresholdat0(2), (slopeat0(1)*slopeat0(2))>0);
    fprintf(fid,'\r\n');
    fclose(fid);
end

disp('(MPFisher) Done.');
return;