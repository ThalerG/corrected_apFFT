# Import necessary libraries
import numpy as np
import matplotlib.pyplot as plt
from corrected_apFFT import corrected_2Win_apFFT

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

# Plotting the result
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