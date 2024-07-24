function [f_correct,y_correct] = Nuttall_2win_apFFT(signal, Fs, thr)
% apFFT Performs the Nuttall double-window all Phase FFT transformation on a 1D array.
% 
%   This function takes a 1D array 'signal' as input and returns its
%   all-phase FFT transformation 'y'. The all-phase FFT aims to minimize
%   phase errors and align all components in phase.
%
%   Inputs:
%       signal - The input 1D array.
%       Fs - The sampling frequency of the signal.
%       thr - The threshold amplitude value for correction. Only peaks with
%             amplitudes greater than thr will be corrected. Default: 0.
%
%   Outputs:
%       f_correct - The corrected frequencies of the peaks.
%       y_correct - The corrected amplitudes and phases of the peaks.
%
%
%   Example usage:
%       [f_correct, y_correct] = Nuttall_2win_apFFT(signal, Fs, thr);

if nargin < 3
    thr = 0; % Default value for thr
end

if mod(length(signal),2)==0
    warning('Signal cannot be of even length. Considering samples 1:(end-1)');
    signal = signal(1:end-1);
end

N = (length(signal)+1)/2;

% Nuttall window
N_win = nuttall(N,'numTerms',4,'Order',3);
N_win = N_win/sum(N_win);

% Convolution of windows (double windowing)
N_win_conv = conv(N_win,N_win);
N_win_conv = N_win_conv/sum(N_win_conv); % Ensures that the sum of the 
                                         % window is 1
% Single windowing of last N samples for FFT
x_win = signal(end-N+1:end).*N_win;
S_FFT = fft(x_win);

% Double-windowing of last 2N-1 samples for apFFT
x_conv = signal.*N_win_conv;
S = [x_conv(N), x_conv(1:N-1) + x_conv(N+1:end)];

S_apFFT = fft(S);

k = 1:floor(N/2); % Index for single-sided FFT/apFFT

 % Amplitude of single-sided FFT/apFFT is doubled for non-zero frequencies

S_apFFT = [S_apFFT(1), 2*S_apFFT(k(2:end))];
S_FFT = [S_FFT(1), 2*S_FFT(k(2:end))];

% Get ApFFT peaks
[pks,locs] = findpeaks(abs(S_apFFT));

locs = locs(pks>thr);

% Frequency and amplitude correction
delta_ph = angle(S_FFT(locs))-angle(S_apFFT(locs));
f_correct = (delta_ph*N/(pi*(N-1))+locs-1)*Fs/N;
A = abs(S_FFT(locs)).^2./abs(S_apFFT(locs));

% Complex representation
ph = angle(S_apFFT(locs));
y_correct = A.*exp(1i.*ph);

end