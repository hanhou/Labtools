% swanalmain.m - matlab program for simulated sine-wave
%             analysis on the simplest lowpass filter:
%
%                 y(n) = x(n)+x(n-1)}

B = [1,1]; % filter feedforward coefficients 
A = 1;     % filter feedback coefficients (none)

N=10;           % number of sinusoidal test frequencies
fs = 1;         % sampling rate in Hz (arbitrary)

fmax = fs/2;    % highest frequency to look at
df = fmax/(N-1);% spacing between frequencies
f = 0:df:fmax;  % sampled frequency axis
dt = 1/fs;      % sampling interval in seconds
tmax = 10;      % number of seconds to run each sine test
t = 0:dt:tmax;  % sampled time axis
ampin = 1;      % test input signal amplitude
phasein = 0;    % test input signal phase

[gains,phases] = swanal(t,f/fs,B,A); % sine-wave analysis

swanalmainplot;    % final plots and comparison to theory