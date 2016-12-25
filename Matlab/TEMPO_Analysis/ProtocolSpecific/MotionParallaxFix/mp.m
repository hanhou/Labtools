function mp(freq, magnitude, d, ex, gain)
% function mp(freq, magnitude, d, ex, gain)
%
% freq (frequency of movement in Hz)
% magnitude (amplitude of movement in m)
% d (depth of patch in deg. of disparity [-]near)
% ex (patch eccentricity x in degrees) USELESS NOW
% gain (pursuit position gain)

%clf;

duration = 2;   % duration (time of trial in s)
vd = 31.9;        % vd (viewing distance in cm)
io = 3.5;       % io (intraocular distance in cm)

ds = -2:1:2
figure;
hold on;

for i = 1:size(ds,2)  % loop for writing out all depths to the trajectories file
%d = ds(i);
    
% dd (depth of patch in cm from fp [+]near)
dd = vd - 180*io*vd/(-d*pi*vd+180*io);  % converts degrees of disparity to depth in cm

% high power gaussian envelope
g_sigma = 55;
g_exponent = 22;
t0 = duration*60/2; 
s = t0/60*g_sigma;

for t=1:duration*1000
    % xt (position of moog in space in cm) = gaussian envelope * sin wave
    xt(t) = exp(-(t*(60/1000)-t0)^g_exponent/(2*s^g_exponent)) * magnitude*100*sin((t/1000)*2*pi*freq);
    % xt(t) = magnitude*100*sin((t/1000)*2*pi*freq);  % w/o gaussian
    xt2(t) = xt(t)*gain;  % Kluge to simulate undergain
    alpha1(t) = atan(vd/xt2(t));
    alpha2(t) = atan((vd-dd)/(xt(t)-((vd-dd)*tan(ex/(180/pi)))));
    alphat(t) = (180/pi)*(alpha1(t)-alpha2(t));  % alphat (retinal angle in degrees)
    if alphat(t) >= 90
        alphat(t) = alphat(t) - 180;
    end
    if alphat(t) < -90
        alphat(t) = alphat(t) + 180;
    end
    if t>1
      dalphat(t) = (alphat(t)-alphat(t-1))/0.001;
    else
      dalphat(t) = 0;
    end
end
plot(xt);
vt = 1000*diff((180/pi)*atan(xt/vd));
plot(vt,'r');
plot(-alphat,'k');
plot(dalphat,'g');
xlabel('Time (ms)');
h = title(sprintf('Freq: %2.2gHz  Magnitude: %gcm\nPatch depth: %gdeg.  Patch eccentricity: %gdeg.',freq,magnitude*100,d,ex));
position = get(h,'Position');
position(2) = position(2)*.75;
set(h,'Position',position);
legend('Monkey position (cm)','Gaze velocity (deg/s)','Retinal angle (deg)','Retinal velocity (deg/s)',3);
max_monkey_speed = max(diff(xt))*1000;
if max_monkey_speed > 50,
    disp('**** Warning!  Monkey moving too fast! ****');
    max_monkey_speed
end
mean_speed_in_degrees_of_retinal_angle = mean(abs(dalphat))

RVt(:,i) = dalphat';  % store for writing out to trajectories file
end  % end loop for writing out to trajectories file

% Write trajectory to a file for import to Origin for figure-making
% PATHOUT = 'Z:\Data\MOOG\Barracuda\Analysis\Trajectories.txt';
% outfile = PATHOUT;
% fid = fopen(outfile, 'w');  % Open text file.
% fprintf(fid, 'Time Xt 180Xt Vt 180Vt EVt 180EVt RVtneg2 180RVtneg2 RVtneg1 180RVtneg1 RVt0 180RVt0 RVt1 180RVt1 RVt2 180RVt2');
% fprintf(fid, '\r\n');
% final_data = [[1:duration*1000]' xt' xt'.*-1 [diff(xt)*1000 0]' [diff(xt)*1000 0]'.*-1 [vt 0]' [vt 0]'.*-1 RVt RVt.*-1];
% for i = 1:duration*1000
%     fprintf(fid,' %+2.4f', final_data(i,:));
%     fprintf(fid,'\r\n');
% end
% fclose(fid);

disp('Done.');