function SpikeRaster(spikes,times,ns)
% Plot spike raster. HH20141116
% spikes: nTrials x nTimes binary matrix

if nargin < 2
    times = 1:size(spikes,2);
    ns = 1:size(spikes,1);
end

[nn,tt] = meshgrid(ns,times);
figure();
plot(tt(logical(spikes')),nn(logical(spikes')),'.k');
