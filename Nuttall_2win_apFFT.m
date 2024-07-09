function [y] = Nuttall_2win_apFFT(signal)
%apFFT Performs the Nuttall double-window all Phase FFT transformation on a 1D array.
%   This function takes a 1D array 'signal' as input and returns its
%   all-phase FFT transformation 'y'. The all-phase FFT aims to minimize
%   phase errors and align all components in phase.

if mod(length(signal),2)==0
    error('Signal cannot be of even length');
end

N = (length(signal)+1)/2;

% Step 1: Nuttall window and N-vector grouping
N_win = nuttall(length(signal),'numTerms',4,'Order',3);
x_win = signal.*N_win;

X_group = NaN([N,N]);

for j = 1:N
    X_group(:,j) = x_win((N-j+1):(end-j+1));
end


% Step 2: Generate the All Phase sequence
% For demonstration, this example assumes a simple phase adjustment.
% In practice, this step might involve more complex phase manipulation
% depending on the specific requirements of the all-phase transformation.
n = length(signal);
phaseFactor = exp(-1i * 2 * pi * (0:n-1) / n);
allPhaseFFT = fftSignal .* phaseFactor;

end