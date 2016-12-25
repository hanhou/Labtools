%-----------------------------------------------------------------------------------------------------------------------
%-- MPBootstrap2.m -- Shuffles (permutes) data and returns p values of original modulations.
%-- Started by JWN, 8/21/06
%-- Last by JWN, 8/21/06
%-----------------------------------------------------------------------------------------------------------------------
function [null_p, fn_p] = MPBootstrap2(data, cond, begin_time, end_time, nulltrue, fntrue);  % 'data' is all data at it's most raw

TEMPO_Defs;
ProtocolDefs;

disp(sprintf('(MPBootstrap2) Started at %s.',datestr(now,14)));

MPdepths = data.moog_params(PATCH_DEPTH,:,MOOG);
uMPdepths = unique(MPdepths);
num_depths = size(uMPdepths,2);
MPtrial_types = data.moog_params(MP_TRIAL_TYPE,:,MOOG);
uMPtrial_types = unique(MPtrial_types);
num_trial_types = length(uMPtrial_types);
MPphase = data.moog_params(MOVEMENT_PHASE,:,MOOG);
uMPphase = unique(MPphase);
num_phase = size(uMPphase,2);
trials = size(MPphase,2);
total_spike_bins = end_time - begin_time;
num_reduced_bins = 39;
bin_width = total_spike_bins/(num_reduced_bins+1);  % ~2000ms/(39+1) = ~50ms;

num_boots = 1000;
booted_null_amps = zeros(1,num_boots);
booted_fn_amps = zeros(1,num_boots);
rand('state',sum(100*clock));  % Reset the uniform generator.

i = cond;
reps = floor(sum(MPtrial_types == i-1)/(num_depths*num_phase));
for boot = 1:num_boots
    for j=[1 2 10]  % Ten MPdepths (which includes null)
        indices = logical((MPtrial_types == i-1) & (MPdepths == uMPdepths(j)));  % Collapsing phases
        raw_spikes = data.spike_data(1,begin_time:end_time,indices);
        shuf = randperm(reps*2);
        hist_data = sum(raw_spikes(1,:,shuf(1:reps)),3);
        [bins, saved_counts0(j,:)] = SpikeBinner(hist_data, 1, bin_width, 0);
        hist_data = sum(raw_spikes(1,:,shuf(reps+1:reps*2)),3);
        [bins, saved_counts180(j,:)] = SpikeBinner(hist_data, 1, bin_width, 0);
    end        
    % First do null (-180 - -0)
    null_counts = saved_counts180(1,:)-saved_counts0(1,:);  % 1 for null
    null_fft = fft(null_counts);
    booted_null_amps(boot) = abs(null_fft(2))/reps;
    % Then do far-near (-)-(2-3)
    a_counts = saved_counts0(10,:)-saved_counts180(2,:);  % 10 and 2 for far and near
    b_counts = saved_counts180(10,:)-saved_counts0(2,:);
    c_counts = b_counts - a_counts; % Aligns phases so that modulation in null that supports far will have the same phase as resulting far-near; see graph paper.  
    fn_fft = fft(c_counts);
    booted_fn_amps(boot) = abs(fn_fft(2))/reps;
end

% Calculate p value
null_p = sum(booted_null_amps > nulltrue)/num_boots;
fn_p = sum(booted_fn_amps > fntrue)/num_boots;

return;