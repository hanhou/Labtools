%-----------------------------------------------------------------------------------------------------------------------
%-- MPBootstrap.m -- Shuffles (permutes) data and returns p value of original PDI.
%-- Then, it bootstraps (resamples with replacement) to generate 95% confidence intervals.
%-- Started by JWN, 9/7/05
%-- Last by JWN, 8/11/06  Changed PDI formula because variances sum.
%-----------------------------------------------------------------------------------------------------------------------
function [p] = MPBootstrap4(data, truePDI, dshift, interpolation_spacing, cleave);  % 'data' can be mean_data or reformatted mis_data

disp(sprintf('(MPBootstrap4) Started at %s.',datestr(now,14)));

num_boots = 1000;
num_depths = size(data,2);
num_reps = size(data,1); % number of data points, which is twice the number of repetitions since we don't care about phase.

mPDI = zeros(1,num_boots);
mGMI = zeros(1,num_boots);
rand('state',sum(100*clock));  % Reset the uniform generator.

% Permutation to decide whether the PDI is sig different from zero.
for i = 1:num_boots
    % Shuffle data
    for k = 1:(num_depths-1)/2  % 'k' is our depths loop incrementer
        rp = randperm(num_reps*2);  
        raw = [data(:,k);data(:,(num_depths+1)-k)];  % Unshuffled column of trials from this depth
        data(:,k) = raw(rp(1:num_reps));  % Pick half and call them near
        data(:,(num_depths+1)-k) = raw(rp(num_reps+1:end));  % Pick other half and call them far
    end

    % Get means and stds of shuffled data
    meanmean_data = mean(data);  % mean of the trials from the means measure
    stdmean_data = std(data);
    % Interpolate the means and stds
    xmeanmean_data = interp1(1:9,meanmean_data,1:interpolation_spacing:9);
    xstdmean_data = interp1(1:9,stdmean_data,1:interpolation_spacing:9);
    % Shift interpolated
    if(dshift > 0)
        sxmeanmean_data = xmeanmean_data(2*cleave+1:length(xmeanmean_data));
        sxstdmean_data = xstdmean_data(2*cleave+1:length(xstdmean_data));
    elseif(dshift <= 0)
        sxmeanmean_data = xmeanmean_data(1:length(xmeanmean_data)-2*cleave);
        sxstdmean_data = xstdmean_data(1:length(xstdmean_data)-2*cleave);
    end
    num_sinterp = length(sxmeanmean_data);
    sxPDI = zeros(1,floor(num_sinterp/2));
    for k = 1:num_sinterp/2  % Using all remaining pairs
        nearm = sxmeanmean_data(k);
        farm = sxmeanmean_data(num_sinterp+1-k);
        nearstd = sxstdmean_data(k);
        farstd = sxstdmean_data(num_sinterp+1-k);
        sxPDI(k) = (farm-nearm)/(abs(farm-nearm)+sqrt((nearstd^2+farstd^2)/2));
    end
    % Save PDI
    mPDI(i) = mean(sxPDI);
end
% Calculate p value
if(truePDI>0)
    p = sum(mPDI>truePDI)/num_boots;
else
    p = sum(mPDI<truePDI)/num_boots;
end

% % Bootstrap to return 95% CIs
% for i = 1:num_boots
%     % Resample data
%     for k = 1:(num_depths-1)/2  % 'k' is our depths loop incrementer
%         crnear = ceil(rand(1,num_reps)*num_reps);
%         near_raw = data(:,k);
%         near_data = near_raw(crnear);
%         crfar = ceil(rand(1,num_reps)*num_reps);
%         far_raw = data(:,(num_depths+1)-k);
%         far_data = far_raw(crfar);
%         % Calculate PDI
%         nearm = mean(near_data);
%         farm = mean(far_data);
%         nearstd = std(near_raw);  % Using fixed stddevs
%         farstd = std(far_raw);
%         PDI(k) = (farm-nearm)/(abs(farm-nearm)+(nearstd+farstd)/2);
%         GMI(k) = (farm-nearm)/(farm+nearm);  % Global modulation index does not use stddev
%     end
%     % Save PDI
%     mPDI(i) = mean(PDI);
%     mGMI(i) = mean(GMI);
% end
% % Kluge to calc trueGMI
% for k = 1:(num_depths-1)/2  % 'k' is our depths loop incrementer
%     near_raw = data(:,k);
%     far_raw = data(:,(num_depths+1)-k);
%     nearm = mean(near_raw);
%     farm = mean(far_raw);
%     nearstd = std(near_raw);  % Using fixed stddevs
%     farstd = std(far_raw);
%     tPDI(k) = (farm-nearm)/(abs(farm-nearm)+(nearstd+farstd)/2);
%     tGMI(k) = (mean(far_raw)-mean(near_raw))/(mean(far_raw)+mean(near_raw));
% end
% trueGMI = mean(tGMI);
% truePDIrecalc = mean(tPDI);
% % Calculate 95% CIs
% smPDI = sort(mPDI);
% truePDI;
% truePDIrecalc;
% cilowD = smPDI(round(num_boots*0.025));
% cihighD = smPDI(round(num_boots*0.975));
% smGMI = sort(mGMI);
% trueGMI;
% cilowM = smGMI(round(num_boots*0.025));
% cihighM = smGMI(round(num_boots*0.975));

return;