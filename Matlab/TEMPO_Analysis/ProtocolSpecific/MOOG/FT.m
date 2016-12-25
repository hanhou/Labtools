function [freq, ampl, phase_angle] = FT(data_x, data_y, FFT_PTS, DC_remove, plot_flag)

DFT_input = data_y - (DC_remove == 1)*mean(data_y);
DFT_mag = abs(fft(DFT_input, FFT_PTS)/length(data_y)); %abs computes the complex modulus when input is complex
DFT = fft(DFT_input, FFT_PTS)/length(data_y);
% added normalization by length of signal (see Matlab Help file on FFT.m, version 2007a or higher) - CRF 4-29-10

phase_angle = angle(DFT);
half_length = floor(length(DFT_mag)/2);
DFT_mag_half = DFT_mag(1:half_length);
phase_angle = phase_angle(1:half_length);
ampl = 2*DFT_mag_half; % multiply by 2 for correct amplitude (to compensate for using half-length?) - CRF 4-29-10

dx = (data_x(length(data_x)) - data_x(1))/(length(data_x)-1);
if dx == 0
    disp('Divide by zero error about to occur in FT.m; place debug point to investigate'); % CRF
    pause;
end
max_freq = 1/(2*dx);
freq = (0:(length(DFT_mag_half)-1)).*(max_freq/length(DFT_mag_half));

if (plot_flag)
    figure;
    subplot(3,1,1);
    bar(data_x, data_y, 'k-');
    title('Input Data');
    
    subplot(3,1,2);
%     plot(freq(1:6), DFT_mag_half(1:6), 'kx-');
    plot(freq, ampl, 'kx-');
    title('Amplitude Spectrum');
    
    subplot(3,1,3);  % added phase plot - CRF 8/17/08
%     plot(freq(1:6), phase_angle(1:6)*180/pi, 'rx-');
    plot(freq, phase_angle*180/pi, 'rx-');
    title('Phase Spectrum');
    xlabel('Frequency (Hz)');
    set(gca,'ytick',[-180 -90 0 90 180]); ylim([-200 200]); 
end

return;