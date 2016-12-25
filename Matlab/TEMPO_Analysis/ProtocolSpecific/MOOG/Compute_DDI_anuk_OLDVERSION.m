% Compute_DDI.m:  Computes a disparity discrimination index.  GCD, 8/11/00
% hdisp should be a ROW vector containing the horizontal disparity of each trial
% resp should be a ROW vector containing the response on each trial
% modified for heading
function [DDI, var_term] = Compute_DDI_anuk(azimuth, elevation, resp)

unique_az = munique(azimuth');
unique_el = munique(elevation');%get unique values of hdisp

%PlotTuningCurve(hdisp, resp, '*','-' , 1, 1)

% %take sqrt of responses to homogenize variance- find negatives first
% negs = ones(1,length(resp));
% negs(find(resp<0)) = -1;
% resp = sqrt(abs(resp));
% resp = resp.*negs;
k=1;
empty = zeros(1,40);
mean_rate = []; stdev = []; resp_ssq = [];
for i=1:length(unique_az)
    for j = 1:length(unique_el)
        matches = (azimuth == unique_az(i)) & (elevation == unique_el(j));
        if sum(matches)>1
            mean_rate(k) = mean( resp(matches) );
            resp_ssq(k) = sum( (resp(matches) - mean_rate(k)).^2 );
        else
            empty(k) = 1;
            mean_rate(k) = 0;
            resp_ssq(k) = 0;
        end
        k = k+1;
    end
end
good = find(empty==0);
max_resp = max(mean_rate(good));
min_resp = min(mean_rate(good));

total_ssq = sum(resp_ssq);
    
var_term = 2*sqrt(total_ssq/(length(resp)-26));
max_min_term = max_resp - min_resp;

DDI = max_min_term/(max_min_term + var_term);

return;