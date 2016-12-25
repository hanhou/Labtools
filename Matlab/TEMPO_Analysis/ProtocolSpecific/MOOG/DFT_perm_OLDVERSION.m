%this takes in 2 sec of data and the original dft ratio of the data. It
%permutes the data and recalculates DFT ratio.
function [p, p_vel, p_acc]= DFT_perm(data, dft_ratio, phase, repeat)
time = 0:.05:1.95;
for i = 1:repeat
    shuf = randperm(length(data));
    shuf_dat = data(shuf);
    [f, amp, resp_phase] = FT(time, shuf_dat, 40, 1, 0);
    f1 = mean(amp(2:4));
    f2 = mean(amp(5:end));
    if f2 == 0
        dft_shuf(i) = 0;
    else
        dft_shuf(i) = f1/f2;
    end
    max_p = resp_phase(find(amp == max(amp(2:4))));
    max_phase = max_p(1);
    dft_vel(i) = -1*dft_shuf(i)*cos(max_phase);
    dft_acc(i) = dft_shuf(i)*sin(max_phase);
end
% p values
p = 2*length(find(dft_shuf > dft_ratio))/repeat;
if -dft_ratio*cos(phase) > 0
    p_vel = 2*length(find(dft_vel > -dft_ratio*cos(phase)))/repeat;
else
    p_vel = 2*length(find(dft_vel < -dft_ratio*cos(phase)))/repeat;
end
if dft_ratio*sin(phase)> 0
    p_acc = 2*length(find(dft_acc > dft_ratio*sin(phase)))/repeat;
else
    p_acc = 2*length(find(dft_acc < dft_ratio*sin(phase)))/repeat;
end
return