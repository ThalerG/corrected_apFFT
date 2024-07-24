# Import necessary libraries
import numpy as np
import matplotlib.pyplot as plt
from corrected_apFFT import Nuttall_2win_apFFT

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