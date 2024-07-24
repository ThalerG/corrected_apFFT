function [f_correct,y_correct] = corrected_2Win_apFFT(signal, Fs, varargin)
% corrected_2Win_apFFT Performs the corrected double-window all Phase FFT transformation on a 1D array.
% 
%   This function takes a 1D array 'signal' as input and returns its
%   all-phase FFT transformation 'y'. The all-phase FFT aims to minimize
%   phase errors and align all components in phase.
%
%   Inputs:
%       signal - The input 1D array.
%       Fs - The sampling frequency of the signal.
%       Name-Value Pairs:
%           'window' - The window function to be applied to the signal before FFT. Default: 'rectwin'.
%           'thr' - The threshold amplitude value for correction. Only peaks with
%                   amplitudes greater than thr will be corrected. Default: 0.
%
%   Outputs:
%       f_correct - The corrected frequencies of the peaks.
%       y_correct - The corrected amplitudes and phases of the peaks.
%
%   Reference:
%       T. Su, M. Yang, T. Jin, and R. C. C. Flesch., "Power harmonic and interharmonic detection method in
%       renewable power based on Nuttall double-window all-phase FFT algorithm", IET Renewable Power 
%       Generation, 12(8), pp 953-961, 2018. https://doi.org/10.1049/iet-rpg.2017.0115
%
%   Example usage:
%       [f_correct, y_correct] = corrected_2Win_apFFT(signal, Fs, 'window', window, 'thr', thr);

% Parse optional inputs
parser = inputParser;
addParameter(parser, 'window', 'rectwin');
addParameter(parser, 'thr', 0);
parse(parser, varargin{:});
window = parser.Results.window;
thr = parser.Results.thr;

if ~ischar(window)
    error('Window must be a string specifying the window function.');
end

if mod(length(signal),2)==0
    warning('Signal for double-windowed apFFT cannot be of even length. Last sample removed.');
    signal = signal(1:end-1);
end

N = (length(signal)+1)/2;

% Select window function
if strcmp(window, 'nuttall4_0')
    win = nuttall(N, 'numTerms', 4, 'Order', 0);
elseif strcmp(window, 'nuttall4_1')
    win = nuttall(N, 'numTerms', 4, 'Order', 1);
elseif strcmp(window, 'nuttall4_3')
    win = nuttall(N, 'numTerms', 4, 'Order', 3);
elseif strcmp(window, 'nuttall4_5')
    win = nuttall(N, 'numTerms', 4, 'Order', 5);
elseif strcmp(window, 'nuttall3_0')
    win = nuttall(N, 'numTerms', 3, 'Order', 0);
elseif strcmp(window, 'nuttall3_1')
    win = nuttall(N, 'numTerms', 3, 'Order', 1);
elseif strcmp(window, 'nuttall3_3')
    win = nuttall(N, 'numTerms', 3, 'Order', 3);
else
    win = generate(sigwin.(window)(N)).';
end

win = win/sum(win); % Ensures that the sum of the 
                    % window is 1

% Convolution of windows (double windowing)
win_conv = conv(win,win);
win_conv = win_conv/sum(win_conv); % Ensures that the sum of the 
                                   % window is 1
% Single windowing of last N samples for FFT
x_win = signal(end-N+1:end).*win;
S_FFT = fft(x_win);

% Double-windowing of last 2N-1 samples for apFFT
x_conv = signal.*win_conv;
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
