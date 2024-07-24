import numpy as np
from nuttall import nuttall
import matplotlib.pyplot as plt
from scipy.fft import fft
from scipy.signal import find_peaks
from scipy.signal.windows import get_window
import warnings
from typing import Tuple

def corrected_2Win_apFFT(signal: np.ndarray, Fs: float, window: str = "boxcar", thr: float = 0) -> Tuple[np.ndarray, np.ndarray]:
    """
    Performs the corrected double-window all-phase FFT transformation on a 1D array.

    Parameters
    ----------
    - signal: 1D array of real numbers representing the input signal.
    - Fs: Sampling frequency of the input signal.
    - window: The window function to be used for the transformation. Default is "boxcar" (rectangular window).
    - thr: Threshold value used to filter out peaks with low amplitude. Default is 0.

    Returns
    -------
    - f_correct: 1D array of real numbers representing the corrected frequencies of the signal components.
    - y_correct: 1D array of complex numbers representing the corrected amplitudes and phases of the signal components.

    This function applies the double-window all-phase FFT transformation to the input signal. The transformation aims to minimize phase errors and align frequency peaks. The function first applies a Nuttall window to the last N samples of the input signal, where N is determined based on the length of the signal. Then, it performs a single-windowed FFT on the windowed samples to obtain the single-sided FFT spectrum. Next, it convolves the window with itself to obtain the double-windowing function, and applies it to the entire input signal. The resulting signal is then transformed using the all-phase FFT. Finally, the function corrects the frequencies and amplitudes of the signal components based on the peaks in the all-phase FFT spectrum.

    Note: The input signal should have an odd length to comply with the requirements of the double-windowed all-phase FFT transformation.
    
    Reference
    ----------
    .. [1] T. Su, M. Yang, T. Jin, and R. C. C. Flesch., "Power harmonic and interharmonic detection method in renewble power based on Nuttall double-window all-phase FFT algorithm", IET Renewable Power Generation, 12(8), pp 953-961, 2018. https://doi.org/10.1049/iet-rpg.2017.0115
    """

    if len(signal) % 2 == 0:
        warnings.warn('Signal for double-windowed apFFT cannot be of even length. Last sample removed.', UserWarning)
        signal = signal[:-1]

    N = (len(signal) + 1) // 2

    
    # Select window function
    if window == "nuttall4_0":
        win = nuttall(N, numTerms=4, Order=0)
    elif window == "nuttall4_1":
        win = nuttall(N, numTerms=4, Order=1)
    elif window == "nuttall4_3":
        win = nuttall(N, numTerms=4, Order=3)
    elif window == "nuttall4_5":
        win = nuttall(N, numTerms=4, Order=5)
    elif window == "nuttall3_0":
        win = nuttall(N, numTerms=3, Order=0)
    elif window == "nuttall3_1":
        win = nuttall(N, numTerms=3, Order=1)
    elif window == "nuttall3_3":
        win = nuttall(N, numTerms=3, Order=3)
    else:
        win = get_window(window, N)

    win /= np.sum(win) # Ensures that the sum of the window is 1

    # Single windowing of last N samples for FFT
    x_win = signal[-N:] * win
    S_FFT = fft(x_win)

    # Convolution of windows (double windowing)
    win_conv = np.convolve(win, win)
    win_conv /= np.sum(win_conv)  # Ensures that the sum of the window is 1

    # Double-windowing of last 2N-1 samples for apFFT
    x_conv = signal * win_conv
    S = np.concatenate(([x_conv[N-1]], x_conv[:N-1] + x_conv[N:]))

    S_apFFT = fft(S)

    # Amplitude of single-sided FFT/apFFT is doubled for non-zero frequencies
    S_apFFT = np.concatenate(([S_apFFT[0]], 2 * S_apFFT[1:N//2]))

    # Amplitude of single-sided FFT/apFFT is doubled for non-zero frequencies
    S_FFT = np.concatenate(([S_FFT[0]], 2 * S_FFT[1:N//2]))

    k = np.arange(0, N//2)  # Index for single-sided FFT/apFFT

    # Correction algorithm only works for peaks
    peaks, _ = find_peaks(abs(S_apFFT[k]))

    # Select only peaks with amplitude greater than the threshold
    peaks = peaks[abs(S_apFFT[peaks]) > thr]

    # Frequency and amplitude correction
    delta_ph = np.angle(S_FFT[peaks]) - np.angle(S_apFFT[peaks])
    f_correct = (delta_ph * N / (np.pi * (N - 1)) + peaks) * Fs / N
    A = np.abs(S_FFT[peaks])**2 / np.abs(S_apFFT[peaks])

    # Complex representation
    ph = np.angle(S_apFFT[peaks])
    y_correct = A * np.exp(1j * ph)
    
    return f_correct, y_correct

if __name__=="__main__":
    # Generate a test signal
    Fs = 1024  # Sampling frequency

    L = 4095 # Length of signal = 2*N-1

    N = (L+1)/2

    f = np.array([50.1, 75.3, 100.4, 125.2, 150.3])  # Frequency of signal components (Hz)
    A = np.array([220, 0.05, 0.1, 25, 7])  # Amplitude of signal components
    ph = np.deg2rad([40, 20, 30, 70, 60])  # Phase of signal components (rad)

    t = np.arange(-N, N+1) / Fs
    x = np.zeros(len(t))

    for k in range(len(f)):
        x += A[k] * np.cos(2 * np.pi * f[k] * t + ph[k])

    # Apply Nuttall_2win_apFFT
    (f_apFFT, y_apFFT) = corrected_2Win_apFFT(x, Fs, window= "nuttall4_3", thr = 0.01)

    plt.figure(figsize=(10, 6))
    ax1 = plt.subplot(2, 1, 1)
    plt.plot(f_apFFT, np.abs(y_apFFT), 'o')
    plt.title('Nuttall double-windowed apFFT of the signal')
    plt.ylabel('Magnitude')
    plt.subplot(2, 1, 2, sharex=ax1)
    plt.plot(f_apFFT, np.angle(y_apFFT, deg=True), 'o')
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('Phase [deg]')
    plt.show()