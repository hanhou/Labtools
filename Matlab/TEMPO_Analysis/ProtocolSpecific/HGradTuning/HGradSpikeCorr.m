function [data] = SpikeCorr(data, SpikeChan, SpikeChan2)

spike_data1 = data.spike_data(SpikeChan, :);
spike_data2 = data.spike_data(SpikeChan2, :);


[correl, lags] = xcorr(spike_data1, spike_data2, 50);
figure
plot(lags, correl);

[val, ind] = max(correl);

choptime = lags(ind);
spikes = find(spike_data1 ~= 0);

%spike_data2(spikes(:) - choptime+1) = spike_data2(spikes(:) - choptime+1)-1;
for i=1:length(spikes)
    if spike_data2(spikes(i) - choptime) ~= 0
        spike_data2(spikes(i) - choptime) = spike_data2(spikes(i) - choptime)-1;
    end
    if spike_data2(spikes(i) - choptime+1) ~= 0
        spike_data2(spikes(i) - choptime+1) = spike_data2(spikes(i) - choptime+1)-1;
    end
    if spike_data2(spikes(i) - choptime-1) ~= 0
        spike_data2(spikes(i) - choptime-1) = spike_data2(spikes(i) - choptime-1)-1;
    end
end

[correl, lags] = xcorr(spike_data1, spike_data2, 50);
hold on
plot(lags, correl, 'r');

data.spike_data(SpikeChan2, :) = spike_data2;