% modified wrapped gaussian w/ 6 params:

% q = [A mu sigma K K-sig DC]

% q(1) = A = amplitude (max firing rate)
% q(2) = mu = position of 1st peak (peaks are always 180 deg apart)
% q(3) = sigma = width of 2nd peak
% q(4) = K = multiplier on 2nd peak amplitude (ratio of peak heights)
% q(5) = K-sig = multiplier on 1st peak width (ratio of peak widths)
% q(6) = DC = DC offset (baseline firing rate)

% mean params for visual:
q_vis = [49.7109 3 2.3821 0.1305 0.6229 5.6826];
% and vestibular:
q_ves = [32.3709 3 2.9991 .17 0.5330 9.2853];

xdata = [0:0.1:360] * pi/180; % (radians)

F_vis = q_vis(1) * ( exp(-2*(1-cos(xdata-q_vis(2)))/(q_vis(5)*q_vis(3))^2) + q_vis(4)*exp(-2*(1-cos(xdata-q_vis(2)-pi))/q_vis(3)^2) ) + q_vis(6);
F_ves = q_ves(1) * ( exp(-2*(1-cos(xdata-q_ves(2)))/(q_ves(5)*q_ves(3))^2) + q_ves(4)*exp(-2*(1-cos(xdata-q_ves(2)-pi))/q_ves(3)^2) ) + q_ves(6);

figure;
plot(xdata, F_vis, 'r', xdata, F_ves, 'k');