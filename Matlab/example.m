clear; close all;

%% Signal:

Fs = 1024; % Sampling freq

L = 4095; % Length of signal = 2*N-1

N = (L+1)/2;

f = [50.1,75.3,100.4,125.2,150.3]; % Frequency of signal components (Hz)
A = [220,0.05,0.1,25,7]; % Amplitude of signal components
ph = deg2rad([40,20,30,70,60]); % Phase of signal components (rad)

% f = [50,80]; % Frequency of signal components (Hz)
% A = [4,2]; % Amplitude of signal components
% ph = [0,-pi/4]; % Phase of signal components (rad)

t=(-N+1:N)/Fs;
x = zeros(1,length(t));

for k = 1:length(f)
    x = x + A(k)*cos(2*pi*f(k)*t+ph(k));
end

% Add Gaussian noise to signal x
noise = 0.01*randn(size(x));
x = x + noise;

%% Corrected double-window apFFT:

% Corrected apFFT
[f_apFFT,y_apFFT] = Nuttall_2win_apFFT(x(1:end-1),Fs);

% Traditional FFT
f_FFT = Fs/L*(0:(L/2));
y_FFT = fft(x)/L;
y_FFT = y_FFT(1:floor(L/2)+1);
y_FFT(2:end-1) = 2*y_FFT(2:end-1);

% Plot comparison
figure;
ax(1) = subplot(2,1,1);
plot(f_FFT,abs(y_FFT));
hold on;
plot(f_apFFT,abs(y_apFFT))
hold off
grid on;
xlabel("Frequency [Hz]")
ylabel("Amplitude")
legend(["FFT","Corrected apFFT"],"Location","best")
ax(2) = subplot(2,1,2);
plot(f_FFT,rad2deg(angle(y_FFT)));
hold on;
plot(f_apFFT,rad2deg(angle(y_apFFT)))
hold off
grid on;
xlabel("Frequency [Hz]")
ylabel("Phase [deg]")
linkaxes(ax,'x')