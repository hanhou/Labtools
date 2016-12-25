function s = noise_canceller(g, freq, harmonics, Fs, mu, printflag)

% ADAPTIVE NOISE CANCELLER - shared by Jing Liu of Newsome lab.
% ------------------------
% Function Call: noice_canceller(g, harmonics, Fs)
% ------------------------------------------------
% g         --> Input Signal 
% freq      --> Fundamental Frequency to be cancelled
% harmonics --> Which harmonics need to be cancelled. Can handle up to 3 harmonics. This should be 
%               a binary vector e.g. [1 0 0] means cancel the 1st harmonic but not the 2nd or 3rd.
% Fs        --> Signal Sampling Frequency
%
% s         --> Output Signal - The function returns the filtered signal.
% printflag --> Whether you want to plot the signals and in time and frequency domains before and 
%               after filtering

harm1 = harmonics(1);
harm2 = harmonics(2);
harm3 = harmonics(3);

nfr   = freq;
% mu    = 0.25;    % Make this smaller if sharper notch is required!

% make the signal zero mean
g = g - mean(g);
%get size of input signal
[M K] = size(g');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adaptive Notch Filtering Implementation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % define initial conditions
    w1 = zeros(1,K); w2 = zeros(1,K);
    w3 = zeros(1,K); w4 = zeros(1,K);
    w5 = zeros(1,K); w6 = zeros(1,K);
    w7 = zeros(1,K); w8 = zeros(1,K);
     e = zeros(1,K);  y = zeros(1,K);

    % define frequency cancellation 
    % waveforms for fundamental freq
    k = 1 : K;
    x1 = cos(2*pi*nfr*k/Fs);
    x2 = sin(2*pi*nfr*k/Fs);

    if harm1 == 1
        x3 = cos(2*pi*nfr*2*k/Fs);
        x4 = sin(2*pi*nfr*2*k/Fs);
    else
        x3(1:K) = 0; x4(1:K) = 0;
        w3(1:K) = 0; w4(1:K) = 0;
    end
    
    if harm2 == 1
        x5 = cos(2*pi*nfr*3*k/Fs);
        x6 = sin(2*pi*nfr*3*k/Fs);
    else
        x5(1:K) = 0; x6(1:K) = 0;
        w5(1:K) = 0; w6(1:K) = 0;
    end
    
    if harm3 == 1
        x7 = cos(2*pi*nfr*4*k/Fs);
        x8 = sin(2*pi*nfr*4*k/Fs);
    else
        x7(1:K) = 0; x8(1:K) = 0;
        w7(1:K) = 0; w8(1:K) = 0;
    end

    % first iteration using initial conditions
    y(1) = w1(1) .* x1(1) + w2(1) .* x2(1) + w3(1) .* x3(1) + w4(1) .* x4(1) + ...
        w5(1) .* x5(1) + w6(1) .* x6(1) + w7(1) .* x7(1) + w8(1) .* x8(1);

    e(1) = g(1) - y(1);

    two_mu = 2 * mu;
        
    for i = 2 : K,

        % update weights
        w1(i) = w1(i-1) + two_mu * (e(i-1) .* x1(i-1));
        w2(i) = w2(i-1) + two_mu * (e(i-1) .* x2(i-1));
        
        if harm1 == 1
            w3(i) = w3(i-1) + two_mu * (e(i-1) .* x3(i-1));
            w4(i) = w4(i-1) + two_mu * (e(i-1) .* x4(i-1));
        end
        
        if harm2 == 1
            w5(i) = w5(i-1) + two_mu * (e(i-1) .* x5(i-1));
            w6(i) = w6(i-1) + two_mu * (e(i-1) .* x6(i-1));
        end
        
        if harm3 == 1
            w7(i) = w7(i-1) + two_mu * (e(i-1) .* x7(i-1));
            w8(i) = w8(i-1) + two_mu * (e(i-1) .* x8(i-1));
        end
        
        % compute adaptive filter output
        y(i) = w1(i) .* x1(i) + w2(i) .* x2(i) + w3(i) .* x3(i) + w4(i) .* x4(i) + ...
            w5(i) .* x5(i) + w6(i) .* x6(i) + w7(i) .* x7(i) + w8(i) .* x8(i);

        % compute noise canceler output
        e(i)  = g(i)  - y(i);
        
   end

   s = e;
      
% end notch filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%run fft on original signal
fo = fftshift(abs(fft(g)));
%get rid of DC 
fo(round(K/2)-5 : round(K/2)+5) = 0;

%run fft on filtered signal
fc = fftshift(abs(fft(s)));
%get rid of DC 
fc(round(K/2)-5 : round(K/2)+5) = 0;

% normalize frequency axis 
F = [-K/2 : K/2-1] / K;
F = Fs.*F;

if printflag
    figure;
    subplot 411; plot(g); title('Original Signal');
    subplot 412; plot(F, fo); title('Original Signal Spectrum'); xlabel('Frequency (Hz)');
    subplot 413; plot(s); title('Filtered Signal');
    subplot 414; plot(F, fc); title('Filtered Signal Spectrum'); xlabel('Frequency (Hz)');
end

