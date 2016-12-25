function y = sine_func(out, x, frq)

% sine_func.m
% calculate y values for a sine function so that you can plot a sine wave
% represented by the output values of suresfit1h_noah.m
% out(1) is the amplitude, or peak displacement of the folded sinusoids
% from the centre.
% x is the neural signal itself (binned spikes)
% out(2) is the phase shift of the sinusoid
% out(3) is the dc offset, or non-zero center amplitude of the sinusoid

y = out(1)*sin(frqf*2*pi*x + (out(2)/180)*pi) + out(3);
return

