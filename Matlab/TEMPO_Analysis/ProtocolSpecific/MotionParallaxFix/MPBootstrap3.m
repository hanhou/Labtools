%-----------------------------------------------------------------------------------------------------------------------
%-- MPBootstrap3.m -- Shuffles time bins (counts) in the final differenced
%-- PSTH.  Good for null_amps or fn_amps.  
%-- Started by JWN, 8/27/06
%-- Last by JWN, 12/17/07  removed reps variable since now I normalize up front
%-----------------------------------------------------------------------------------------------------------------------
function [p] = MPBootstrap3(counts, ffttrue);  % counts come from final differenced PSTH

disp(sprintf('(MPBootstrap3) Started at %s.',datestr(now,14)));

num_boots = 1000;
num_bins = length(counts);
fft_amps = zeros(1,num_boots);
rand('state',sum(100*clock));  % Reset the uniform generator.

% Permutation to decide whether the fft amplitude is sig different from zero.
for i = 1:num_boots
    % Shuffle data
    shuffled_counts = randsample(counts,num_bins);
    fft_ = fft(shuffled_counts);
    % Save fft amplitude (1st harmonic only)
    fft_amps(i) = abs(fft_(2));
end
% Calculate p value
p = sum(fft_amps > ffttrue)/num_boots;

return;