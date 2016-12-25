resp = load('Z:\LabTools\Matlab\TEMPO_Analysis\ProtocolSpecific\SpeedTuning\SpeedTuningRawRates_Tania.dat');
speed = [0 0.5 1 2 4 8 16 32]; 

select_resp = resp(:,1:8);

peak_resps = max(select_resp')';

norm_resp = zeros(size(select_resp));
for i=1:length(peak_resps)
    norm_resp(i,:) = select_resp(i,:)/peak_resps(i);
end

pop_mean = mean(norm_resp);
pop_se = std(norm_resp)/sqrt(length(peak_resps));

[speed' pop_mean' pop_se']

figure;
errorbar(speed, pop_mean, pop_se);
ylim([0 1]);
xlabel('Speed (deg/sec)');
ylabel('Normalized Response');
title('MT Population Speed Tuning Curve (n=499)');