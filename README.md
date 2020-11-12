# Zaf-Julia

Zafar's Audio Functions in Julia for audio signal analysis (UNDER CONSTRUCTION).
- [`zaf.py`](#zafpy): Julia module with the audio functions.
- [`examples.ipynb`](#examplesipynb): Jupyter module with some examples.
- [`audio_file.wav`](#audio_filewav): audio file used for the examples.

## zaf.py

This Julia module implements a number of functions for audio signal analysis. Simply copy the file `zaf.jl` in your working directory and you are good to go. Make sure to have the following packages installed (`Pkg.add("name_of_the_package")`):
- [WAV](https://juliapackages.com/p/wav): Julia package to read and write the WAV audio file format.
- [FFTW](https://juliapackages.com/p/fftw): Julia bindings to the [FFTW](http://www.fftw.org/) library for fast Fourier transforms (FFTs), as well as functionality useful for signal processing.
- [Plots](https://docs.juliaplots.org/latest/): powerful convenience for visualization in Julia.

Functions:
- `stft` - [Short-time Fourier transform (STFT)](#short-time-fourier-transform-stft)
- `istft` - [Inverse STFT](#inverse-short-time-fourier-transform-stft)
- `cqtkernel` - [Constant-Q transform (CQT) kernel](#constant-q-transform-cqt-kernel)
- `cqtspectrogram` - [CQT spectrogram using a CQT kernel](#constant-q-transform-cqt-spectrogram-using-a-cqt-kernel)
- `cqtchromagram` - [CQT chromagram using a CQT kernel](#constant-q-transform-cqt-chromagram-using-a-cqt-kernel)
- `mfcc` - [Mel frequency cepstrum coefficients (MFCCs)](#mel-frequency-cepstrum-coefficients-mfccs)
- `dct` - [Discrete cosine transform (DCT) using the fast Fourier transform (FFT)](#discrete-cosine-transform-dct-using-the-fast-fourier-transform-fft)
- `dst` - [Discrete sine transform (DST) using the FFT](#discrete-sine-transform-dst-using-the-fast-fourier-transform-fft)
- `mdct` - [Modified discrete cosine transform (MDCT) using the FFT](#modified-discrete-cosine-transform-mdct-using-the-fast-fourier-transform-fft)
- `imdct` - [Inverse MDCT using the FFT](#inverse-modified-discrete-cosine-transform-mdct-using-the-fast-fourier-transform-fft)

Other:
- `sigplot` - Plot an audio signal in seconds
- `specshow` - Display an audio spectrogram in dB, seconds, and Hz
- `cqtspecshow` - Display a CQT audio spectrogram in dB, seconds, and Hz
- `cqtchromshow` - Display a CQT audio chromagram in seconds


### Short-time Fourier transform (STFT)

```
audio_stft = zaf.stft(audio_signal, window_function, step_length)
    
Inputs:
    audio_signal: audio signal (number_samples,)
    window_function: window function (window_length,)
    step_length: step length in samples
Output:
    audio_stft: audio STFT (window_length, number_frames)
```

#### Example: compute and display the spectrogram from an audio file

```
# Load the modules
include("./zaf.jl")
using .zaf
using WAV
using Statistics
using Plots

# Read the audio signal (normalized) with its sampling frequency in Hz, and average it over its channels
audio_signal, sampling_frequency = wavread("audio_file.wav");
audio_signal = mean(audio_signal, dims=2);

# Set the window duration in seconds (audio is stationary around 40 milliseconds)
window_duration = 0.04;

# Derive the window length in samples (use powers of 2 for faster FFT and constant overlap-add (COLA))
window_length = nextpow(2, ceil(Int, window_duration*sampling_frequency));

# Compute the window function (periodic Hamming window for COLA)
window_function = zaf.hamming(window_length, "periodic");

# Set the step length in samples (half of the window length for COLA)
step_length = convert(Int, window_length/2);

# Compute the STFT
audio_stft = zaf.stft(audio_signal, window_function, step_length)

# Derive the magnitude spectrogram (without the DC component and the mirrored frequencies)
audio_spectrogram = abs.(audio_stft[2:convert(Int, window_length/2)+1, :])

# Magnitude spectrogram (without the DC component and the mirrored frequencies)
audio_spectrogram = abs.(audio_stft[2:convert(Int, window_length/2)+1,:]);

# Display the spectrogram in dB, seconds, and Hz
xtick_step = 1
ytick_step = 1000
zaf.specshow(audio_spectrogram, length(audio_signal), sampling_frequency, xtick_step, ytick_step)
heatmap!(title = "Spectrogram (dB)", size = (990, 600))
```


<img src="images/stft.png" width="1000">


### Inverse short-time Fourier transform (STFT)

```
audio_signal = zaf.istft(audio_stft, window_function, step_length)

Inputs:
    audio_stft: audio STFT (window_length, number_frames)
    window_function: window function (window_length,)
    step_length: step length in samples
Output:
    audio_signal: audio signal (number_samples,)
```

#### Example: estimate the center and the sides from a stereo audio file

```
# Import the modules
```


## examples.ipynb

This Jupyter notebook shows some examples for the different functions of the Julia module `zaf`.

See [Jupyter notebook viewer](https://nbviewer.jupyter.org/github/zafarrafii/Zaf-Julia/blob/master/examples.ipynb).


## audio_file.wav

23 second audio excerpt from the song *Que Pena Tanto Faz* performed by *Tamy*.


# Author

- Zafar Rafii
- zafarrafii@gmail.com
- http://zafarrafii.com/
- [CV](http://zafarrafii.com/Zafar%20Rafii%20-%20C.V..pdf)
- [GitHub](https://github.com/zafarrafii)
- [LinkedIn](https://www.linkedin.com/in/zafarrafii/)
- [Google Scholar](https://scholar.google.com/citations?user=8wbS2EsAAAAJ&hl=en)
