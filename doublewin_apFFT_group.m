function [f_correct,y_correct] = doublewin_apFFT_group(signal, Fs)
%apFFT Performs the Nuttall double-window all Phase FFT transformation on a 1D array.
%   This function takes a 1D array 'signal' as input and returns its
%   all-phase FFT transformation 'y'. The all-phase FFT aims to minimize
%   phase errors and align all components in phase.

if mod(length(signal),2)==0
    error('Signal cannot be of even length');
end

plot_graphs = true;

N = (length(signal)+1)/2;
tol = 1e-6;

% Nuttall window
N_win = nuttall(N,'numTerms',4,'Order',3);
N_win = N_win/sum(N_win);

% Nuttall window
NAll_win = nuttall(length(signal),'numTerms',4,'Order',3);
NAll_win = NAll_win/sum(NAll_win);

% Single windowing of last N samples for FFT
x_win = signal(end-N+1:end).*N_win;
S_FFT = fft(x_win);
S_FFT(abs(S_FFT)<tol) = 0;

X_group = NaN([N,N]);
x_signal_win = signal.*NAll_win;

for j = 1:N
    X_group(:,j) = x_signal_win((N-j+1):(end-j+1));
    X_group(:,j) = circshift(X_group(:,j),-j+1);
end

S = mean(X_group,2)';

S_apFFT = fft(S);

k = 1:floor(N/2); % Index for single-sided FFT/apFFT

 % Amplitude of single-sided FFT/apFFT is doubled for non-zero frequencies

S_apFFT = [S_apFFT(1), 2*S_apFFT(k(2:end))];
S_FFT = [S_FFT(1), 2*S_FFT(k(2:end))];

% Frequency and amplitude correction
delta_ph = angle(S_FFT(k))-angle(S_apFFT(k));
f_correct = (delta_ph*N/(pi*(N-1))+k-1)*Fs/N;
A = abs(S_FFT(k)).^2./abs(S_apFFT(k));

% Complex representation
ph = angle(S_apFFT(k));
y_correct = A.*exp(1i.*ph);

if plot_graphs
    f_FFT = Fs/N*(0:(N/2));
    
    figure;
    ax(1) = subplot(3,1,1);
    plot(f_FFT(1:end-1),abs(S_FFT(k)),'LineWidth',1)
    hold on
    plot(f_FFT(1:end-1),abs(S_apFFT(k)),'LineWidth',1)
    legend("FFT","apFFT"); ylabel("Amplitude");
    hold off;
    ax(2) = subplot(3,1,2);
    plot(f_FFT(1:end-1),angle(S_FFT(k)),'LineWidth',1)
    hold on
    plot(f_FFT(1:end-1),angle(S_apFFT(k)),'LineWidth',1)
    legend("FFT","apFFT"); ylabel("Phase [rad]");
    hold off;
    ax(3) = subplot(3,1,3);
    plot(f_correct,abs(A),'LineWidth',1)
    legend("Corrected apFFT"); ylabel("Amplitude");
    linkaxes(ax,'x')
    xlabel("Frequency [Hz]")
end

end