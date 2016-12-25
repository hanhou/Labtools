%-----------------------------------------------------------------------------------------------------------------------
%-- Vergence_ChoiceProb.m -- Uses ROC analysis to compute a h_verg choice probability for each different stimulus level
%--	GCD, 2/14/01
%-----------------------------------------------------------------------------------------------------------------------
function Vergence_ChoiceProb(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

TEMPO_Defs;		%needed for defines like IN_T1_WIN_CD
ProtocolDefs;	%needed for all protocol specific functions - contains keywords - BJP 1/4/01

Pref_HDisp = data.one_time_params(PREFERRED_HDISP);

%get the column of values of horiz. disparities in the dots_params matrix
h_disp = data.dots_params(DOTS_HDISP,:,PATCH1);
unique_hdisp = munique(h_disp');

%get the binocular correlations
binoc_corr = data.dots_params(DOTS_BIN_CORR, :, PATCH1);
unique_bin_corr = munique(binoc_corr');

%get signed binocular correlations
sign = (h_disp == Pref_HDisp)*2 - 1;	%=1 if preferred disparity, -1 if null disparity
signed_bin_corr = binoc_corr .* sign;
unique_signed_bin_corr = munique(signed_bin_corr');

%now, get the firing rates for all the trials 
spike_rates = data.spike_rates(SpikeChan, :);

%compute the raw and calibrated eye position data for this run
[h_verg, v_verg, h_conj, v_conj, calib_h_verg, calib_v_verg, calib_h_conj, calib_v_conj] = ...
    DDiscrim_GetEyeData(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

%get indices of any NULL conditions (for measuring spontaneous activity
null_trials = logical( (binoc_corr == data.one_time_params(NULL_VALUE)) );

%now, determine the choice that was made for each trial, PREFERRED or NULL
%by definition, a preferred choice will be made to Target1 and a null choice to Target 2
%thus, look for the events IN_T1_WIN_CD and IN_T2_WIN_CD.  GCD, 5/30/2000
num_trials = length(binoc_corr);
PREFERRED = 1;
NULL = 2;
for i=1:num_trials
    temp = data.event_data(1,:,i);
    events = temp(temp>0);  % all non-zero entries
    if (sum(events == IN_T1_WIN_CD) > 0)
        choice(i) = PREFERRED;
    elseif (sum(events == IN_T2_WIN_CD) > 0)
        choice(i) = NULL;
    else
        disp('Neither T1 or T2 chosen.  This should not happen!.  File must be bogus.');
    end        
end

%prepare data for ANCOVA analysis
binoc_corr_ranks = binoc_corr;
for i=1:length(unique_bin_corr)
    select = (binoc_corr == unique_bin_corr(i));
    binoc_corr_ranks(select) = i;
end
hdisp_ranks = h_disp;
for i=1:length(unique_hdisp)
    select = (h_disp == unique_hdisp(i));
    hdisp_ranks(select) = i;
end

%now, Z-score the spike rates for each bin_corr and disparity condition
Z_Spikes = spike_rates;
for i=1:length(unique_bin_corr)
    for j=1:length(unique_hdisp)
        select = (binoc_corr == unique_bin_corr(i)) & (h_disp == unique_hdisp(j));
        z_dist = spike_rates(select);
        z_dist = (z_dist - mean(z_dist))/std(z_dist);
        Z_Spikes(select) = z_dist;
    end
end

%[binoc_corr_ranks' hdisp_ranks' choice' spike_rates' h_verg']

%run the ANCOVA analysis tool
aoctool(calib_h_verg', Z_Spikes', choice);

Z_Spikes_pref = Z_Spikes;
Z_Spikes_pref(choice == 2) = NaN;
Z_Spikes_null = Z_Spikes;
Z_Spikes_null(choice == 1) = NaN;
Mdata = [calib_h_verg' Z_Spikes_pref' Z_Spikes_null' Z_Spikes' choice'];

i = size(PATH,2) - 1;
while PATH(i) ~='\'	%Analysis directory is one branch below Raw Data Dir
    i = i - 1;
end   
PATHOUT = [PATH(1:i) 'Analysis\NeuroPsychoCurves\'];
i = size(FILE,2) - 1;
while FILE(i) ~='.'
    i = i - 1;
end
FILEOUT = [FILE(1:i) 'resp_verg_ANCOVA'];
fileid = [PATHOUT FILEOUT];
fwriteid = eval(['fopen(fileid, ''w'')']);

for i=1:length(hdisp_ranks)
    fprintf(fwriteid, '%4d %4d %4d %7.3f %8.5f\n', binoc_corr_ranks(i), hdisp_ranks(i), choice(i), spike_rates(i), h_verg(i) );
end
fclose(fwriteid);

return;