function [freq, ampl, phase_angle] = FT(data_x, data_y, FFT_PTS, DC_remove, plot_flag)

DFT_input = data_y - (DC_remove == 1)*mean(data_y);
DFT_mag = abs(fft(DFT_input, FFT_PTS)); %abs computes the complex modulus when input is complex
DFT = fft(DFT_input, FFT_PTS);

phase_angle = angle(DFT);
half_length = floor(length(DFT_mag)/2);
DFT_mag_half = DFT_mag(1:half_length);
phase_angle = phase_angle(1:half_length);
ampl = DFT_mag_half;

dx = (data_x(length(data_x)) - data_x(1))/(length(data_x)-1);
if dx == 0
    disp('Divide by zero error about to occur in FT.m; place debug point.');
    pause;
end
max_freq = 1/(2*dx);
freq = (0:(length(DFT_mag_half)-1)).*(max_freq/length(DFT_mag_half));

if (plot_flag)
    figure;
    subplot(2,1,1);
    bar(data_x, data_y, 'k-');
    title('Input Data');
    
    subplot(2,1,2);
    plot(freq, DFT_mag_half, 'kx-');
    title('Amplitude Spectrum');
    xlabel('Frequency (Hz)');
end

return;