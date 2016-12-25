function [gains,phases] = swanal(t,f,B,A)
% SWANAL - Perform sine-wave analysis on filter B(z)/A(z)
  
ampin = 1;      % input signal amplitude
phasein = 0;    % input signal phase

N = length(f);       % number of test frequencies
gains = zeros(1,N);  % pre-allocate amp-response array
phases = zeros(1,N); % pre-allocate phase-response array
if length(A)==1
  ntransient=length(B)-1; % no. samples to steady state
else
  error('Need to set transient response duration here');
end

for k=1:length(f)    % loop over analysis frequencies
  s = ampin*cos(2*pi*f(k)*t+phasein); % test sinusoid
  y = filter(B,A,s); % run it through the filter
  yss = y(ntransient+1:length(y)); % chop off transient
  % measure output amplitude as max (SHOULD INTERPOLATE):
  [ampout,peakloc] = max(abs(yss)); % ampl. peak & index
  gains(k) = ampout/ampin;  % amplitude response
  if ampout < eps  % eps returns "machine epsilon"
    phaseout=0;    % phase is arbitrary at zero gain
  else
    sphase = 2*pi*f(k)*(peakloc+ntransient-1);
    % compute phase by inverting sinusoid (BAD METHOD):
    phaseout = acos(yss(peakloc)/ampout) - sphase;
    phaseout = mod2pi(phaseout); % reduce to [-pi,pi)
  end
  phases(k) = phaseout-phasein;
  swanalplot; % signal plotting script
end