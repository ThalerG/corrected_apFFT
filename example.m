clear; close all;

%% Signal:

Fs = 1024; % Sampling freq
N = 2048; % Number of samples

f = [50.1,75.3,100.4,125.2,150.3]; % Frequency of signal components (Hz)
A = [220,0.05,0.1,25,7]; % Amplitude of signal components
ph = deg2rad([40,20,30,70,60]); % Phase of signal components (rad)

% f = [50,80]; % Frequency of signal components (Hz)
% A = [4,2]; % Amplitude of signal components
% ph = [0,-pi/4]; % Phase of signal components (rad)

t=(0:N-1)/Fs;
x = zeros(1,length(t));

for k = 1:length(f)
    x = x + A(k)*cos(2*pi*f(k)*t+ph(k));
end

%% Corrected double-window apFFT:

% Corrected apFFT
[f_apFFT,y_apFFT] = Nuttall_2win_apFFT(x(1:end-1),Fs);
% [f_apFFT2,y_apFFT2] = doublewin_apFFT_group(x(1:end-1),Fs);

% Traditional FFT
f_FFT = Fs/N*(0:(N/2));
y_FFT = fft(x)/N;
y_FFT = y_FFT(1:N/2+1);
y_FFT(2:end-1) = 2*y_FFT(2:end-1);

figure;
plot(f_FFT,abs(y_FFT));
hold on;
plot(f_apFFT,abs(y_apFFT))
hold off
xlabel("Frequency [Hz]")
ylabel("Amplitude")