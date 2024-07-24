import numpy as np
from windows import nuttall
import matplotlib.pyplot as plt
from scipy.fft import fft

def Nuttall_2win_apFFT(signal, Fs):
    """
    Performs the Nuttall double-window all Phase FFT transformation on a 1D array.
    This function takes a 1D array 'signal' as input and returns its
    all-phase FFT transformation 'y'. The all-phase FFT aims to minimize
    phase errors and align all components in phase.
    """

    CONV = False # Convolution of windows (double windowing) or by groups

    if len(signal) % 2 == 0:
        raise ValueError('Signal cannot be of even length')

    N = (len(signal) + 1) // 2

    # Nuttall window
    N_win = nuttall(N,4,3)
    # N_win = np.ones((N))
    N_win /= np.sum(N_win)    

    ################### BY CONVOLUTION ############################

    # Convolution of windows (double windowing)
    N_win_conv = np.convolve(N_win, N_win)
    N_win_conv /= np.sum(N_win_conv)  # Ensures that the sum of the window is 1

    # Double-windowing of last 2N-1 samples for apFFT
    x_conv = signal * N_win_conv
    S = np.concatenate(([x_conv[N-1]], x_conv[:N-1] + x_conv[N:]))

    S_apFFT_conv = fft(S)

    # Amplitude of single-sided FFT/apFFT is doubled for non-zero frequencies
    S_apFFT_conv = np.concatenate(([S_apFFT_conv[0]], 2 * S_apFFT_conv[1:N//2]))

    ################### BY GROUPS ############################
    
    # Nuttall window for the whole signal
    NAll_win = nuttall(len(signal), 4, 3)
    NAll_win /= np.sum(NAll_win)

    X_group = np.empty((N, N))
    X_group[:] = np.NaN

    # Step 1.1: First windowing
    x_signal_win = signal*NAll_win

    for j in range(N):
        # Step 1.2: Grouping in N-length vectors
        x_temp = x_signal_win[(N - j - 1):(len(signal) - j)]

        # Step 2: Circular shifting
        X_group[j, :] = np.roll(x_temp, -j)

    # Step 3: Transpose matrix
    X_group = X_group.T

    # Step 4.1: Second windowing
    X_group_w = X_group*N_win

    # Step 4.2: Sum of the N-length vectors
    S_gr = np.sum(X_group_w,axis = 1)

    # Step 4.3: FFT
    S_apFFT_gr = fft(S_gr)

    # Amplitude of single-sided FFT/apFFT is doubled for non-zero frequencies
    S_apFFT_gr = np.concatenate(([S_apFFT_gr[0]], 2 * S_apFFT_gr[1:N//2]))

    diff = np.abs(S_gr - S)

    time = np.arange(len(diff)) / Fs

    plt.figure(figsize=(10, 6))

    # Plot difference and elementwise division
    plt.subplot(3, 1, 1)
    plt.plot(time, diff)
    plt.title('Difference')
    plt.xlabel('Time (s)')
    plt.ylabel('Difference between methods')

    plt.subplot(3, 1, 2)
    plt.plot(time, S_gr / S)
    plt.xlabel('Time (s)')
    plt.ylabel('Ratio between methods')

    # Plot signals S_gr and S over time
    ax1 = plt.subplot(3, 1, 3)
    ax2 = ax1.twinx()
    ax1.plot(time, S_gr, label='S_gr', color='blue')
    ax2.plot(time, S, label='S', color='red')
    ax1.set_xlabel('Time (s)')
    ax1.set_ylabel('Amplitude (antes)', color='blue')
    ax2.set_ylabel('Amplitude (depois)', color='red')
    ax1.legend(loc='upper left')
    ax2.legend(loc='upper right')

    plt.show()

    # Single windowing of last N samples for FFT
    x_win = signal[-N:] * N_win
    S_FFT = fft(x_win)

    # Amplitude of single-sided FFT/apFFT is doubled for non-zero frequencies
    S_FFT = np.concatenate(([S_FFT[0]], 2 * S_FFT[1:N//2]))

    k = np.arange(0, N//2)  # Index for single-sided FFT/apFFT

    f = np.linspace(0, Fs/2, len(S_FFT))

    if CONV:
        S_apFFT = S_apFFT_conv
    else:
        S_apFFT = S_apFFT_gr

    # Frequency and amplitude correction
    delta_ph = np.angle(S_FFT) - np.angle(S_apFFT)
    f_correct = (delta_ph * N / (np.pi * (N - 1)) + k) * Fs / N
    A = np.abs(S_FFT)**2 / np.abs(S_apFFT)

    # Complex representation
    ph = np.angle(S_apFFT)
    y_correct = A * np.exp(1j * ph)


    plt.figure(figsize=(10, 6))

    # Plot magnitude
    plt.subplot(2, 1, 1)
    plt.plot(f, np.abs(S_apFFT_conv), label='S_apFFT (antes)')
    # plt.plot(f, np.abs(S_apFFT_gr), label='S_apFFT (depois)')
    # plt.plot(f, np.abs(S_FFT), label='S_FFT')
    plt.title('Nuttall 2win apFFT of the signal')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Magnitude')
    plt.legend()

    # Plot phase
    plt.subplot(2, 1, 2)
    plt.plot(f, np.angle(S_apFFT_conv, deg = True), label='ApFFT (antes)')
    # plt.plot(f, np.angle(S_apFFT_gr, deg = True), label='ApFFT (depois)')
    # plt.plot(f, np.angle(S_FFT, deg = True), label='FFT')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Phase (deg)')
    plt.legend()

    plt.show()
    
    return f_correct, y_correct

if __name__=="__main__":
    # Generate a test signal
    Fs = 1024  # Sampling frequency
    N = 2047  # Number of samples, odd to comply with Nuttall_2win_apFFT requirements

    f = np.array([50.1, 75.3, 100.4, 125.2, 150.3])  # Frequency of signal components (Hz)
    A = np.array([220, 0.05, 0.1, 25, 7])  # Amplitude of signal components
    ph = np.deg2rad([40, 20, 30, 70, 60])  # Phase of signal components (rad)

    t = np.arange(N) / Fs
    x = np.zeros(len(t))

    for k in range(len(f)):
        x += A[k] * np.cos(2 * np.pi * f[k] * t + ph[k])

    # Apply Nuttall_2win_apFFT
    (f_apFFT, y_apFFT) = Nuttall_2win_apFFT(x, Fs)

    # Plotting the result
    plt.figure(figsize=(10, 6))
    plt.plot(f_apFFT, np.abs(y_apFFT), label='Nuttall 2win apFFT')
    plt.title('Nuttall 2win apFFT of the signal')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Magnitude')
    plt.legend()
    plt.show()