# Advanced Phase FFT Analysis with Nuttall Windowing

This repository contains MATLAB and Python implementations for correcting peak location in frequency domain using all phase Fast Fourier Transform (FFT) and double windowing. The primary focus is on minimizing phase errors and aligning all components in phase for a given signal.

This repository is based on the following research paper:

Su, T., Yang, M., Jin, T., & Flesch, R. C. C. (2018). Power harmonic and interharmonic detection method in renewable power based on Nuttall double-window all-phase FFT algorithm. IET Renewable Power Generation, 12(8), 953-961. DOI: [10.1049/iet-rpg.2017.0115]( https://doi.org/10.1049/iet-rpg.2017.0115)

## Structure

The repository is organized into several directories, each containing relevant files for MATLAB and Python implementations, as well as test examples and utility functions.

### Directory Structure

- `Matlab/`: Contains MATLAB scripts and functions for the double windowed apFFT analysis and nth order Nuttall window.
  - [`corrected_2Win_apFFT.m`](Matlab/corrected_2Win_apFFT.m): Main function for performing the corrected double-window all Phase FFT transformation.
  - [`nuttall.m`](Matlab/nuttall.m): Function for generating nth order Nuttall window.
  - [`example.m`](Matlab/example.m): Example script demonstrating the use of `corrected_2Win_apFFT` in Matlab using Nuttall window.
- `Python/`: Contains Python scripts for a similar analysis using NumPy and Matplotlib.
  - [`corrected_apFFT.py`](Python/corrected_apFFT.py): Python implementation of the corrected double windowed apFFT.
  - [`windows.py`](Python/windows.py): Function for generating nth order Nuttall window.
  - [`example.py`](Python/example.py): Example script demonstrating the use of `corrected_2Win_apFFT` in Python using Nuttall window.

### Key Files

- MATLAB Implementation: [`corrected_2Win_apFFT.m`](Matlab/corrected_2Win_apFFT.m)
- Python Implementation: [`corrected_apFFT.py`](Python/corrected_apFFT.py)
- Example Usage: [`example.py`](Python/example.py) (Python), [`example.m`](Matlab/example.m) (MATLAB)

## Getting Started

To use the MATLAB functions, ensure you have MATLAB installed and add the `Matlab/` directory to your MATLAB path. For Python, ensure you have Python installed along with NumPy and Matplotlib libraries.

### MATLAB Example

```matlab
% Load your signal
signal = load('your_signal.mat');
Fs = 1024; % Sampling frequency
win = 'nuttall4_3' % 4 terms, 3rd order Nuttall window
thr = 0.001; % Threshold amplitude for correction

% Perform Nuttall windowed apFFT
[f_correct, y_correct] = corrected_2Win_apFFT(signal, Fs, 'window', win, 'thr', thr);
```

### Python Example

```python
# Import the function
from corrected_apFFT import corrected_2Win_apFFT

# Generate or load your signal
Fs = 1024  # Sampling frequency
signal = your_signal_generator_or_loader()

# Perform Nuttall windowed apFFT
f_apFFT, y_apFFT = corrected_2Win_apFFT(signal, Fs,  window= "nuttall4_3", thr=0.001)

# Plot the result
# Refer to `example.py` for a complete example
```

## Contributing
Contributions to improve the implementations or extend the functionality are welcome. Please feel free to fork the repository, make your changes, and submit a pull request.

## License
This project is open-source and available under the MIT License.
