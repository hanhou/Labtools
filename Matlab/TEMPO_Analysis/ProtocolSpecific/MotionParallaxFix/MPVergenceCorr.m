%-----------------------------------------------------------------------------------------------------------------------
%-- MPVergenceCorr.m -- Correlates mean difference in eye position (L-R) with simulated depth in the MP, BD, RM, and C conditions.
% Could be easily expnaded to include analysis of EO and HO conditions.
%-- Started by JWN, 11/16/08
%-- Last by JWN, 11/16/08
%-----------------------------------------------------------------------------------------------------------------------
function MPVergenceCorr(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE);

ver = '1.0';
TEMPO_Defs;
Path_Defs;
symbols = {'bo' 'rs' 'gd' 'kv' 'm<' 'c>' 'bv' 'rv'};
line_types2 = {'b--' 'r--' 'g--' 'k--' 'g.-' 'b.-' 'r-.' 'k.'};
line_types4 = {'b-' 'r-' 'g-' 'k-' 'm-' 'c-' 'y-' 'b-'};
line_types5 = {'bo-' 'rs-' 'gd-' 'kv-' 'm<-' 'c>-' 'yo-' 'bs-'};
NULL_VALUE = -9999;

disp(sprintf('(MPVergenceCorr v%s) Started at %s.',ver,datestr(now,14)));

% Get the trial type, depth values, and movement phase for each condition in the condition_list[]
MPdepths = data.moog_params(PATCH_DEPTH,:,MOOG);
uMPdepths = unique(MPdepths);
num_depths = size(uMPdepths,2);
MPtrial_types = data.moog_params(MP_TRIAL_TYPE,:,MOOG);
uMPtrial_types = unique(MPtrial_types);  % Conditions present
%Place breakouts here.  This is what SelectiveAnalysis could be all about!
%if(isempty(find(uMPtrial_types==###))) return;  end;
if(isempty(find(uMPtrial_types==1))) disp('(MPVergenceCorr) Breakout: No BD');  return;  end;  % BREAKOUT ENABLED!

num_trial_types = length(uMPtrial_types);
MPphase = data.moog_params(MOVEMENT_PHASE,:,MOOG);
uMPphase = unique(MPphase);
num_phase = size(uMPphase,2);
trials = size(MPphase,2);

no2=0;
area='MT';
[monkid, cellid, runstr]=strread(FILE,'m%dc%dr%s.htb');
% Eye check based on file name (set up for Barracuda and Ovid 050107)
switch monkid
    case 9, % Barracuda
        if cellid==155, no2 = 1; end
        if cellid==156, no2 = 1; end
    case 15, % Ovid
        if cellid<35, no2 = 1; end
        if cellid>46, no2 = 1; end
    otherwise
        disp('(MPVergenceCorr) WARNING: Unknown monkid');
end
if no2==1, return; end
% In data.eye_data, Channels 1,2,3&4 are eye (x&y), 5&6 are Moog (x&y).
% Only analyze stimulus time 214:614 (2s long).
eye_xl = data.eye_data(1,215:614,:);
eye_xr = data.eye_data(3,215:614,:);
eye_diff = squeeze(mean(eye_xl-eye_xr));  % Difference between of the two eyes (a measure of vergence)

for i = 1:4
    % Now get indices for BD non-null trials
    indices = logical((MPtrial_types == i-1) & (MPdepths ~= -9999));
    % Get r and p
    [r(i),p(i)] = corr(eye_diff(indices),MPdepths(indices)');
end

PATHOUT = 'Z:\Data\MOOG\Ovid\Analysis\';
filenames = {'VergCorr'};
for i = 1:1
    outfile = cell2mat(strcat(PATHOUT,area,'_',filenames(i),'.txt'));
    headerflag = 0;
    if (exist(outfile) == 0) % File does not yet exist, so print a header
        headerflag = 1;
    end
    fid = fopen(outfile, 'a');  % Open text file.
    if (headerflag)
        fprintf(fid, 'FILE ');
        fprintf(fid, 'MPr BDr RMr Cr MPp BDp RMp Cp');
        fprintf(fid, '\r\n');
    end
    fprintf(fid,'%10s', strtok(FILE,'.'));
    fprintf(fid,' %+2.5f', r, p);
    fprintf(fid,'\r\n');
    fclose(fid);
end

disp('(MPVergenceCorr) Done.');
return;