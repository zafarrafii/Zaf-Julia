"""
This Julia module implements a number of functions for audio signal analysis.

Functions:
    stft - Short-time Fourier transform (STFT)
    istft - inverse STFT
    cqtkernel - Constant-Q transform (CQT) kernel
    cqtspectrogram - CQT spectrogram using a CQT kernel
    cqtchromagram - CQT chromagram using a CQT kernel
    mfcc - Mel frequency cepstrum coefficients (MFCCs)
    dct - Discrete cosine transform (DCT) using the fast Fourier transform (FFT)
    dst - Discrete sine transform (DST) using the FFT
    mdct - Modified discrete cosine transform (MDCT) using the FFT
    imdct - Inverse MDCT using the FFT

Other:
    wavread - Read a WAVE file (using Scipy)
    wavwrite - Write a WAVE file (using Scipy)
    sigplot - Plot an audio signal in seconds
    specshow - Display an audio spectrogram in dB, seconds, and Hz
    cqtspecshow - Display a CQT audio spectrogram in dB, seconds, and Hz
    cqtchromshow - Display a CQT audio chromagram in seconds

Author:
    Zafar Rafii
    zafarrafii@gmail.com
    http://zafarrafii.com
    https://github.com/zafarrafii
    https://www.linkedin.com/in/zafarrafii/
    11/04/20
"""

module zaf

export stft, istft, cqtkernel, cqtspectrogram, cqtchromagram, mfcc, dct, dst,
mdct, imdct

"""
    audio_stft = z.stft(audio_signal, window_function, step_length);

Compute the short-time Fourier transform (STFT)

# Arguments:
- `audio_signal::Float`: the audio signal [number_samples, 1]
- `window_function::Float`: the window function [window_length, 1]
- `step_length::Integer`: the step length in samples
- `audio_stft::Complex`: the audio STFT [window_length, number_frames]

# Example: Compute the spectrogram of an audio file
```
# Audio signal averaged over its channels and sample rate in Hz
Pkg.add("WAV")
using WAV
audio_signal, sample_rate = wavread("audio_file.wav");
audio_signal = mean(audio_signal, 2);

# Window duration in seconds (audio is stationary around 40 milliseconds)
window_duration = 0.04;

# Window length in samples (power of 2 for fast FFT and constant overlap-add (COLA))
window_length = nextpow2(ceil(Int64, window_duration*sample_rate));

# Window function (periodic Hamming window for COLA)
include("z.jl")
window_function = z.hamming(window_length, "periodic");

# Step length in samples (half the window length for COLA)
step_length = convert(Int64, window_length/2);

# Magnitude spectrogram (without the DC component and the mirrored frequencies)
audio_stft = z.stft(audio_signal, window_function, step_length);
audio_spectrogram = abs.(audio_stft[2:convert(Int64, window_length/2)+1,:]);

# Spectrogram displayed in dB, s, and kHz
Pkg.add("Plots")
using Plots
plotly()
x_labels = [string(round(i*step_length/sample_rate, 2)) for i = 1:size(audio_spectrogram, 2)];
y_labels = [string(round(i*sample_rate/window_length/1000, 2)) for i = 1:size(audio_spectrogram, 1)];
heatmap(x_labels, y_labels, 20*log10.(audio_spectrogram))
```
"""
function stft(audio_signal, window_function, step_length)

    # Get the number of samples and the window length in samples
    number_samples = length(audio_signal)
    window_length = length(window_function)

    # Number of time number_frames
    number_times = ceil(Int64, (window_length-step_length+number_samples)/step_length)

    # Zero-padding at the start and end to center the windows
    audio_signal = [zeros(window_length-step_length,1); audio_signal;
    zeros(number_times*step_length-number_samples,1)]

    # Initialize the STFT
    audio_stft = zeros(window_length, number_times)

    # Loop over the time frames
    for time_index = 1:number_times

        # Window the signal
        sample_index = step_length*(time_index-1)
        audio_stft[:, time_index] = audio_signal[1+sample_index:window_length+sample_index].*window_function

    end

    # Fourier transform of the frames
    audio_stft = fft(audio_stft, 1)

end

"""
    audio_istft = z.istft(audio_signal, window_function, step_length);

Compute the inverse short-time Fourier transform (STFT)

# Arguments:
- `audio_stft::Complex`: the audio STFT [window_length, number_frames]
- `window_function::Float`: the window function [window_length, 1]
- `step_length::Integer`: the step length in samples
- `audio_signal::Float`: the audio signal [number_samples, 1]

# Example: Estimate the center and sides signals of a stereo audio file
```
# Stereo audio signal and sample rate in Hz
Pkg.add("WAV")
using WAV
audio_signal, sample_rate = wavread("audio_file.wav");

# Parameters for the STFT
include("z.jl")
window_duration = 0.04;
window_length = nextpow2(ceil(Int64, window_duration*sample_rate));
window_function = z.hamming(window_length,"periodic");
step_length = convert(Int64, window_length/2);

# STFT of the left and right channels
audio_stft1 = z.stft(audio_signal[:,1], window_function, step_length);
audio_stft2 = z.stft(audio_signal[:,2], window_function, step_length);

# Magnitude spectrogram (with DC component) of the left and right channels
audio_spectrogram1 = abs.(audio_stft1[1:convert(Int64, window_length/2)+1, :]);
audio_spectrogram2 = abs.(audio_stft2[1:convert(Int64, window_length/2)+1, :]);

# Time-frequency masks of the left and right channels for the center signal
center_mask1 = min.(audio_spectrogram1, audio_spectrogram2)./audio_spectrogram1;
center_mask2 = min.(audio_spectrogram1, audio_spectrogram2)./audio_spectrogram2;

# STFT of the left and right channels for the center signal (with extension to mirrored frequencies)
center_stft1 = [center_mask1; center_mask1[convert(Int64, window_length/2):-1:2,:]].*audio_stft1;
center_stft2 = [center_mask2; center_mask2[convert(Int64, window_length/2):-1:2,:]].*audio_stft2;

# Synthesized signals of the left and right channels for the center signal
center_signal1 = z.istft(center_stft1, window_function, step_length);
center_signal2 = z.istft(center_stft2, window_function, step_length);

# Final stereo center and sides signals
center_signal = [center_signal1, center_signal2];
center_signal = center_signal[1:size(audio_signal, 1), :];
sides_signal = audio_signal-center_signal;

# Synthesized center and side signals
wavwrite(center_signal, "center_signal.wav", Fs=sample_rate);
wavwrite(sides_signal, "sides_signal.wav", Fs=sample_rate);

# Spectrogram displayed in dB, s, and kHz
Pkg.add("Plots")
using Plots
plotly()
time_signal = (1:size(audio_signal, 1))/sample_rate;
audio_plot = plot(time_signal, audio_signal, xlabel="Time (s)", title="Original Signal");
center_plot = plot(time_signal, center_signal, xlabel="Time (s)", title="Center Signal");
sides_plot = plot(time_signal, sides_signal, xlabel="Time (s)", title="Sides Signal");
plot(audio_plot, center_plot, sides_plot, layout=(3,1), legend=false)
```
"""
function istft(audio_stft, window_function, step_length)

    # Window length in samples and number of time frames
    window_length, number_times = size(audio_stft);

    # Number of samples for the signal
    number_samples = (number_times-1)*step_length+window_length;

    # Initialize the signal
    audio_signal = zeros(number_samples, 1);

    # Inverse Fourier transform of the frames and real part to ensure real values
    audio_stft = real(ifft(audio_stft, 1));

    # Loop over the time frames
    for time_index = 1:number_times

        # Constant overlap-add (if proper window and step)
        sample_index = step_length*(time_index-1);
        audio_signal[1+sample_index:window_length+sample_index] = audio_signal[1+sample_index:window_length+sample_index]
        + audio_stft[:,time_index];

    end

    # Remove the zero-padding at the start and end
    audio_signal = audio_signal[window_length-step_length+1:number_samples-(window_length-step_length)];

    # Un-apply window (just in case)
    audio_signal = audio_signal/sum(window_function[1:step_length:window_length]);

end

"""
    cqt_kernel = z.cqtkernel(sample_rate, frequency_resolution, minimum_frequency, maximum_frequency);

Compute the constant-Q transform (CQT) kernel

# Arguments:
- `sample_rate::Float` the sample rate in Hz
- `frequency_resolution::Integer` the frequency resolution in number of frequency channels per semitone
- `minimum_frequency::Float`: the minimum frequency in Hz
- `maximum_frequency::Float`: the maximum frequency in Hz
- `cqt_kernel::Complex`: the CQT kernel [number_frequencies, fft_length]

# Example: Compute and display the CQT kernel
```
# CQT kernel parameters
sample_rate = 44100;
frequency_resolution = 2;
minimum_frequency = 55;
maximum_frequency = sample_rate/2;

# CQT kernel
include("z.jl")
cqt_kernel = z.cqtkernel(sample_rate, frequency_resolution, minimum_frequency, maximum_frequency);

# Magnitude CQT kernel displayed
Pkg.add("Plots")
using Plots
plotly()
heatmap(abs.(cqt_kernel))
```
"""
function cqtkernel(sample_rate, frequency_resolution, minimum_frequency, maximum_frequency)

    # Number of frequency channels per octave
    octave_resolution = 12*frequency_resolution;

    # Constant ratio of frequency to resolution (= fk/(fk+1-fk))
    quality_factor = 1/(2^(1/octave_resolution)-1);

    # Number of frequency channels for the CQT
    number_frequencies = round(Int64, octave_resolution*log2(maximum_frequency/minimum_frequency));

    # Window length for the FFT (= window length of the minimum frequency = longest window)
    fft_length = nextpow2(ceil(Int64, quality_factor*sample_rate/minimum_frequency));

    # Initialize the kernel
    cqt_kernel = zeros(Complex64, number_frequencies, fft_length);

    # Loop over the frequency channels
    for frequency_index = 1:number_frequencies

        # Frequency value (in Hz)
        frequency_value = minimum_frequency*2^((frequency_index-1)/octave_resolution);

        # Window length (nearest odd value because the complex exponential will have an odd length, in samples)
        window_length = 2*round(Int64, quality_factor*sample_rate/frequency_value/2)+1;

        # Temporal kernel (without zero-padding, odd and symmetric)
        temporal_kernel = hamming(window_length, "symmetric").*
        exp.(2*pi*im*quality_factor*(-(window_length-1)/2:(window_length-1)/2)/window_length)/window_length;

        # Pre and post zero-padding to center FFTs
        temporal_kernel = [zeros(1, convert(Int64, (fft_length-window_length+1)/2)),
        temporal_kernel', zeros(1, convert(Int64, (fft_length-window_length-1)/2))];

        # Spectral kernel (mostly real because temporal kernel almost symmetric)
        # (Note that Julia's fft equals the complex conjugate of Matlab's fft!)
        spectral_kernel = fft(temporal_kernel);

        # Save the spectral kernels
        cqt_kernel[frequency_index, :] = spectral_kernel;

    end

    # Energy threshold for making the kernel sparse
    energy_threshold = 0.01;

    # Make the CQT kernel sparser
    cqt_kernel[abs.(cqt_kernel).<energy_threshold] = 0;

    # Make the CQT kernel sparse
    cqt_kernel = sparse(cqt_kernel);

    # From Parseval's theorem
    cqt_kernel = conj(cqt_kernel)/fft_length;

end

"""
    audio_spectrogram = z.cqtspectrogram(audio_signal,sample_rate,time_resolution,cqt_kernel);

Compute constant-Q transform (CQT) spectrogram using a kernel

# Arguments:
- `audio_signal::Float`: the audio signal [number_samples, 1]
- `sample_rate::Float`: the sample rate in Hz
- `time_resolution::Float`: the time resolution in number of time frames per second
- `cqt_kernel::Complex`: the CQT kernel [number_frequencies, fft_length]
- `audio_spectrogram::Float`: the audio spectrogram in magnitude [number_frequencies, number_times]

# Example: Compute and display the CQT spectrogram
```
# Audio file averaged over the channels and sample rate in Hz
Pkg.add("WAV")
using WAV
audio_signal, sample_rate = wavread("audio_file.wav");
audio_signal = mean(audio_signal, 2);

# CQT kernel
frequency_resolution = 2;
minimum_frequency = 55;
maximum_frequency = 3520;
include("z.jl")
cqt_kernel = z.cqtkernel(sample_rate, frequency_resolution, minimum_frequency, maximum_frequency);

# CQT spectrogram
time_resolution = 25;
audio_spectrogram = z.cqtspectrogram(audio_signal, sample_rate, time_resolution, cqt_kernel);

# CQT spectrogram displayed in dB, s, and Hz
Pkg.add("Plots")
using Plots
plotly()
x_labels = [string(round(i/time_resolution, 2)) for i = 1:size(audio_spectrogram, 2)];
y_labels = [string(round(55*2^((i-1)/(12*frequency_resolution)), 2)) for i = 1:size(audio_spectrogram, 1)];
heatmap(x_labels, y_labels, 20*log10.(audio_spectrogram))
```
"""
function cqtspectrogram(audio_signal, sample_rate, time_resolution, cqt_kernel)

    # Number of time samples per time frame
    step_length = round(Int64, sample_rate/time_resolution);

    # Number of time frames
    number_times = floor(Int64, length(audio_signal)/step_length);

    # Number of frequency channels and FFT length
    number_frequencies, fft_length = size(cqt_kernel);

    # Zero-padding to center the CQT
    audio_signal = [zeros(ceil(Int64, (fft_length-step_length)/2),1); audio_signal;
    zeros(floor(Int64, (fft_length-step_length)/2),1)];

    # Initialize the spectrogram
    audio_spectrogram = zeros(number_frequencies, number_times);

    # Loop over the time frames
    for time_index = 1:number_times

        # Magnitude CQT using the kernel
        sample_index = step_length*(time_index-1);
        audio_spectrogram[:, time_index] = abs.(cqt_kernel * fft(audio_signal[sample_index+1:sample_index+fft_length]));

    end

    return audio_spectrogram

end

"""
    audio_chromagram = z.cqtchromagram(audio_signal, sample_rate, time_resolution, frequency_resolution, cqt_kernel);

Compute the constant-Q transform (CQT) chromagram using a kernel

# Arguments:
- `audio_signal::Float`: the audio signal [number_samples, 1]
- `sample_rate::Float`: the sample rate in Hz
- `time_resolution::Float`: the time resolution in number of time frames per second
- `frequency_resolution::Integer`: the frequency resolution in number of frequency channels per semitones
- `cqt_kernel::Complex`: the CQT kernel [number_frequencies, fft_length]
- `audio_chromagram::Complex`: the audio chromagram [number_chromas, number_times]

# Example: Compute and display the CQT chromagram
```
# Audio file averaged over the channels and sample rate in Hz
Pkg.add("WAV")
using WAV
audio_signal, sample_rate = wavread("audio_file.wav");
audio_signal = mean(audio_signal, 2);

# CQT kernel
frequency_resolution = 2;
minimum_frequency = 55;
maximum_frequency = 3520;
include("z.jl")
cqt_kernel = z.cqtkernel(sample_rate, frequency_resolution, minimum_frequency, maximum_frequency);

# CQT chromagram
time_resolution = 25;
audio_chromagram = z.cqtchromagram(audio_signal, sample_rate, time_resolution, frequency_resolution, cqt_kernel);

# CQT chromagram displayed in dB, s, and chromas
Pkg.add("Plots")
using Plots
plotly()
x_labels = [string(round(i/time_resolution, 2)) for i = 1:size(audio_chromagram, 2)];
y_labels = [string(i) for i = 1:size(audio_chromagram, 1)];
heatmap(x_labels, y_labels, 20*log10.(audio_chromagram))
```
"""
function cqtchromagram(audio_signal, sample_rate, time_resolution, frequency_resolution, cqt_kernel)

    # CQT spectrogram
    audio_spectrogram = z.cqtspectrogram(audio_signal, sample_rate, time_resolution, cqt_kernel);

    # Number of frequency channels and time frames
    number_frequencies, number_times = size(audio_spectrogram);

    # Number of chroma bins
    number_chromas = 12*frequency_resolution;

    # Initialize the chromagram
    audio_chromagram = zeros(number_chromas, number_times);

    # Loop over the chroma bins
    for chroma_index = 1:number_chromas

        # Sum the energy of the frequency channels for every chroma
        audio_chromagram[chroma_index, :] = sum(audio_spectrogram[chroma_index:number_chromas:number_frequencies, :], 1);

    end

    return audio_chromagram

end

"""
audio_mfcc = z.mfcc(audio_signal, sample_rate, number_filters, number_coefficients);

    Compute the mel frequency cepstrum coefficients (MFFCs)

# Arguments:
- `audio_signal::Float`: the audio signal [number_samples, 1]
- `sample_rate::Float`: the sample rate in Hz
- `number_filters::Integer`: the number of filters
- `number_coefficients::Integer`: the number of coefficients (without the 0th coefficient)
- `audio_mfcc::Float`: the audio MFCCs [number_times, number_coefficients]

# Example: Compute and display the MFCCs, delta MFCCs, and delta-detla MFCCs
```
# Audio signal averaged over its channels and sample rate in Hz
Pkg.add("WAV")
using WAV
audio_signal, sample_rate = wavread("audio_file.wav");
audio_signal = mean(audio_signal, 2);

#  MFCCs for a given number of filters and coefficients
number_filters = 40;
number_coefficients = 20;
include("z.jl")
audio_mfcc = z.mfcc(audio_signal, sample_rate, number_filters, number_coefficients);

# Delta and delta-delta MFCCs
audio_deltamfcc = diff(audio_mfcc, 2);
audio_deltadeltamfcc = diff(audio_deltamfcc, 2);

# MFCCs, delta MFCCs, and delta-delta MFCCs displayed in s
Pkg.add("Plots")
using Plots
plotly()
step_length = convert(Int64, nextpow2(ceil(Int64, 0.04*sample_rate))/2);
time_signal = round.((1:size(audio_mfcc, 2))*step_length/sample_rate, 2);
mfcc_plot = plot(time_signal, audio_mfcc', xlabel="Time (s)", title="MFCCs");
deltamfcc_plot = plot(time_signal[2:end], audio_deltamfcc', xlabel="Time (s)", title="Delta MFCCs");
deltadeltamfcc_plot = plot(time_signal[3:end], audio_deltadeltamfcc', xlabel="Time (s)", title="Delta-delta MFCCs");
plot(mfcc_plot, deltamfcc_plot, deltadeltamfcc_plot, layout=(3,1), legend=false)
```
"""
function mfcc(audio_signal, sample_rate, number_filters, number_coefficients)

    # Window duration in seconds, length in samples, and function, and step length in samples
    window_duration = 0.04;
    window_length = nextpow2(ceil(Int64, window_duration*sample_rate));
    window_function = z.hamming(window_length, "periodic");
    step_length = convert(Int64, window_length/2);

    # Magnitude spectrogram (without the DC component and the mirrored frequencies)
    audio_stft = z.stft(audio_signal, window_function, step_length);
    audio_spectrogram = abs.(audio_stft[2:convert(Int64, window_length/2)+1, :]);

    # Minimum and maximum mel frequencies
    mininum_melfrequency = 2595*log10(1+(sample_rate/window_length)/700);
    maximum_melfrequency = 2595*log10(1+(sample_rate/2)/700);

    # Indices of the overlapping filters (linearly spaced in the mel scale and logarithmically spaced in the linear scale)
    filter_width = 2*(maximum_melfrequency-mininum_melfrequency)/(number_filters+1);
    filter_indices = mininum_melfrequency:filter_width/2:maximum_melfrequency;
    filter_indices = round.(Int64, 700*(10.^(filter_indices/2595)-1)*window_length/sample_rate);

    # Initialize the filter bank
    filter_bank = zeros(number_filters, convert(Int64, window_length/2));

    # Loop over the filters
    for filter_index = 1:number_filters

        # Left and right sides of the triangular overlapping filters (linspace more accurate than triang or bartlett!)
        filter_bank[filter_index,filter_indices[filter_index]:filter_indices[filter_index+1]] =
        linspace(0, 1, filter_indices[filter_index+1]-filter_indices[filter_index]+1);
        filter_bank[filter_index,filter_indices[filter_index+1]:filter_indices[filter_index+2]] =
        linspace(1, 0, filter_indices[filter_index+2]-filter_indices[filter_index+1]+1);

    end

    # Discrete cosine transform of the log of the magnitude spectrogram mapped onto the mel scale using the filter bank
    audio_mfcc = dct(log.(filter_bank*audio_spectrogram+eps()), 1);

    # The first coefficients (without the 0th) represent the MFCCs
    audio_mfcc = audio_mfcc[2:number_coefficients+1, :];

end

"""
audio_dct = z.dct(audio_signal, dct_type);

    Compute the discrete cosine transform (DCT) using the fast Fourier transform (FFT)

# Arguments:
- `audio_signal::Float`: the audio signal [number_samples, number_frames]
- dct_type::Integer`: the DCT type (1, 2, 3, or 4)
- audio_dct::Float`: the audio DCT [number_frequencies, number_frames]

# Example: Compute the 4 different DCTs and compare them to Julia's DCTs
```
# Audio signal averaged over its channels and sample rate in Hz
Pkg.add("WAV")
using WAV
audio_signal, sample_rate = wavread("audio_file.wav");
audio_signal = mean(audio_signal, 2);

# Audio signal for a given window length, and one frame
window_length = 1024;
audio_signal = audio_signal[1:window_length, :];

# DCT-I, II, III, and IV
include("z.jl")
audio_dct1 = z.dct(audio_signal, 1);
audio_dct2 = z.dct(audio_signal, 2);
audio_dct3 = z.dct(audio_signal, 3);
audio_dct4 = z.dct(audio_signal, 4);

# Julia's DCT-II (Julia does not have a DCT-I, III, and IV!)
julia_dct2 = dct(audio_signal, 1);

# DCT-I, II, III, and IV, Julia's version, and their errors displayed
Pkg.add("Plots")
using Plots
plotly()
dct1_plot = plot(audio_dct1, title="DCT-I");
dct2_plot = plot(audio_dct2, title="DCT-II");
dct3_plot = plot(audio_dct3, title="DCT-III");
dct4_plot = plot(audio_dct4, title="DCT-IV");
jdct1_plot = plot(zeros(window_length, 1))
jdct2_plot = plot(audio_dct2, title="Julia's DCT-II");
jdct3_plot = plot(zeros(window_length, 1));
jdct4_plot = plot(zeros(window_length, 1));
zjdct1_plot = plot(zeros(window_length, 1));
zjdct2_plot = plot(audio_dct2-julia_dct2, title="Error");
zjdct3_plot = plot(zeros(window_length, 1));
zjdct4_plot = plot(zeros(window_length, 1));
zeros_plot = plot(zeros(window_length, 1));
plot(dct1_plot, jdct1_plot, zjdct1_plot, dct2_plot, jdct2_plot, zjdct2_plot,
dct3_plot, jdct3_plot, zjdct3_plot, dct4_plot, jdct4_plot, zjdct4_plot, layout=(4,3), legend=false)
```
"""
function dct(audio_signal, dct_type)

    if dct_type == 1

        # Number of samples per frame
        window_length = size(audio_signal, 1);

        # Pre-processing to make the DCT-I matrix orthogonal (concatenate to avoid the input to change!)
        audio_signal = [audio_signal[1:1, :]*sqrt(2); audio_signal[2:window_length-1, :];
        audio_signal[window_length:window_length, :]*sqrt(2)];

        # Compute the DCT-I using the FFT
        audio_dct = [audio_signal; audio_signal[window_length-1:-1:2, :]];
        audio_dct = fft(audio_dct, 1);
        audio_dct = real(audio_dct[1:window_length, :])/2;

        # Post-processing to make the DCT-I matrix orthogonal
        audio_dct[[1, window_length], :] = audio_dct[[1, window_length], :]/sqrt(2);
        audio_dct = audio_dct*sqrt(2/(window_length-1));

    elseif dct_type == 2

        # Number of samples and frames
        window_length, number_frames = size(audio_signal);

        # Compute the DCT-II using the FFT
        audio_dct = zeros(4*window_length, number_frames);
        audio_dct[2:2:2*window_length, :] = audio_signal;
        audio_dct[2*window_length+2:2:4*window_length, :] = audio_signal[window_length:-1:1, :];
        audio_dct = fft(audio_dct, 1);
        audio_dct = real(audio_dct[1:window_length,:])/2;

        # Post-processing to make the DCT-II matrix orthogonal
        audio_dct[1, :] = audio_dct[1, :]/sqrt(2);
        audio_dct = audio_dct*sqrt(2/window_length);

    elseif dct_type == 3

        # Number of samples and frames
        window_length, number_frames = size(audio_signal);

        # Pre-processing to make the DCT-III matrix orthogonal (concatenate to avoid the input to change!)
        audio_signal = [audio_signal[1:1, :]*sqrt(2); audio_signal[2:window_length, :]];

        # Compute the DCT-III using the FFT
        audio_dct = zeros(4*window_length, number_frames);
        audio_dct[1:window_length, :] = audio_signal;
        audio_dct[window_length+2:2*window_length+1, :] = -audio_signal[window_length:-1:1, :];
        audio_dct[2*window_length+2:3*window_length, :] = -audio_signal[2:window_length, :];
        audio_dct[3*window_length+2:4*window_length, :] = audio_signal[window_length:-1:2, :];
        audio_dct = fft(audio_dct, 1);
        audio_dct = real(audio_dct[2:2:2*window_length, :])/4;

        # Post-processing to make the DCT-III matrix orthogonal
        audio_dct = audio_dct*sqrt(2/window_length);

    elseif dct_type == 4

        # Number of samples and frames
        window_length, number_frames = size(audio_signal);

        # Compute the DCT-IV using the FFT
        audio_dct = zeros(8*window_length, number_frames);
        audio_dct[2:2:2*window_length, :] = audio_signal;
        audio_dct[2*window_length+2:2:4*window_length, :] = -audio_signal[window_length:-1:1, :];
        audio_dct[4*window_length+2:2:6*window_length, :] = -audio_signal;
        audio_dct[6*window_length+2:2:8*window_length, :] = audio_signal[window_length:-1:1, :];
        audio_dct = fft(audio_dct, 1);
        audio_dct = real(audio_dct[2:2:2*window_length, :])/4;

        # Post-processing to make the DCT-IV matrix orthogonal
        audio_dct = audio_dct*sqrt(2/window_length);

    end

end

"""
audio_dst = z.dst(audio_signal, dst_type);

    Compute the discrete sine transform (DST) using the fast Fourier transform (FFT)

# Arguments:
- `audio_signal::Float`: the audio signal [number_samples, number_frames]
- `dst_type::Integer`: the DST type (1, 2, 3, or 4)
- `audio_dst::Float`: the audio DST [number_frequencies, number_frames]

# Example: Compute the 4 different DSTs and compare them to their respective inverses
```
# Audio signal averaged over its channels and sample rate in Hz
Pkg.add("WAV")
using WAV
audio_signal, sample_rate = wavread("audio_file.wav");
audio_signal = mean(audio_signal, 2);

# Audio signal for a given window length, and one frame
window_length = 1024;
audio_signal = audio_signal[1:window_length, :];

# DST-I, II, III, and IV
include("z.jl")
audio_dst1 = z.dst(audio_signal, 1);
audio_dst2 = z.dst(audio_signal, 2);
audio_dst3 = z.dst(audio_signal, 3);
audio_dst4 = z.dst(audio_signal, 4);

# Respective inverses, i.e., DST-I, II, III, and IV
audio_idst1 = z.dst(audio_dst1, 1);
audio_idst2 = z.dst(audio_dst2, 3);
audio_idst3 = z.dst(audio_dst3, 2);
audio_idst4 = z.dst(audio_dst4, 4);

# DST-I, II, III, and IV, respective inverses, and errors displayed
Pkg.add("Plots")
using Plots
plotly()
dst1_plot = plot(audio_dst1, title="DST-I");
dst2_plot = plot(audio_dst2, title="DST-II");
dst3_plot = plot(audio_dst3, title="DST-III");
dst4_plot = plot(audio_dst4, title="DST-IV");
idst1_plot = plot(audio_idst1, title="Inverse DST-I = DST-I");
idst2_plot = plot(audio_idst2, title="Inverse DST-II = DST-III");
idst3_plot = plot(audio_idst3, title="Inverse DST-III = DST-II");
idst4_plot = plot(audio_idst4, title="Inverse DST-IV = DST-IV");
ddst1_plot = plot(audio_signal-audio_idst1, title="Error");
ddst2_plot = plot(audio_signal-audio_idst2, title="Error");
ddst3_plot = plot(audio_signal-audio_idst3, title="Error");
ddst4_plot = plot(audio_signal-audio_idst4, title="Error");
plot(dst1_plot, idst1_plot, ddst1_plot, dst2_plot, idst2_plot, ddst2_plot,
dst3_plot, idst3_plot, ddst3_plot, dst4_plot, idst4_plot, ddst4_plot, layout=(4,3), legend=false)
```
"""
function dst(audio_signal, dst_type)

    if dst_type==1

        # Number of samples per frame
        window_length, number_frames = size(audio_signal);

        # Compute the DST-I using the FFT
        audio_dst = [zeros(1, number_frames); audio_signal;
        zeros(1, number_frames); -audio_signal[window_length:-1:1, :]];
        audio_dst = fft(audio_dst, 1);
        audio_dst = -imag(audio_dst[2:window_length+1, :])/2;

        # Post-processing to make the DST-I matrix orthogonal
        audio_dst = audio_dst*sqrt(2/(window_length+1));

    elseif dst_type==2

        # Number of samples per frame
        window_length, number_frames = size(audio_signal);

        # Compute the DST-II using the FFT
        audio_dst = zeros(4*window_length, number_frames);
        audio_dst[2:2:2*window_length, :] = audio_signal;
        audio_dst[2*window_length+2:2:4*window_length, :] = -audio_signal[window_length:-1:1, :];
        audio_dst = fft(audio_dst, 1);
        audio_dst = -imag(audio_dst[2:window_length+1, :])/2;

        # Post-processing to make the DST-II matrix orthogonal
        audio_dst[window_length, :] = audio_dst[window_length, :]/sqrt(2);
        audio_dst = audio_dst*sqrt(2/window_length);

    elseif dst_type==3

        # Number of samples per frame
        window_length, number_frames = size(audio_signal);

        # Pre-processing to make the DST-III matrix orthogonal (concatenate to avoid the input to change!)
        audio_signal = [audio_signal[1:window_length-1, :];
        audio_signal[window_length:window_length, :]*sqrt(2)];

        # Compute the DST-III using the FFT
        audio_dst = zeros(4*window_length, number_frames);
        audio_dst[2:window_length+1, :] = audio_signal;
        audio_dst[window_length+2:2*window_length, :] = audio_signal[window_length-1:-1:1, :];
        audio_dst[2*window_length+2:3*window_length+1, :] = -audio_signal;
        audio_dst[3*window_length+2:4*window_length, :] = -audio_signal[window_length-1:-1:1, :];
        audio_dst = fft(audio_dst, 1);
        audio_dst = -imag(audio_dst[2:2:2*window_length, :])/4;

        # Post-processing to make the DST-III matrix orthogonal
        audio_dst = audio_dst*sqrt(2/window_length);

    elseif dst_type==4

        # Number of samples per frame
        window_length, number_frames = size(audio_signal);

        # Compute the DST-IV using the FFT
        audio_dst = zeros(8*window_length, number_frames);
        audio_dst[2:2:2*window_length, :] = audio_signal;
        audio_dst[2*window_length+2:2:4*window_length, :] = audio_signal[window_length:-1:1, :];
        audio_dst[4*window_length+2:2:6*window_length, :] = -audio_signal;
        audio_dst[6*window_length+2:2:8*window_length, :] = -audio_signal[window_length:-1:1, :];
        audio_dst = fft(audio_dst, 1);
        audio_dst = -imag(audio_dst[2:2:2*window_length, :])/4;

        # Post-processing to make the DST-IV matrix orthogonal
        audio_dst = audio_dst*sqrt(2/window_length);

    end

end

"""
audio_mdct = z.mdct(audio_signal,window_function);

    Compute the modified discrete cosine transform (MDCT) using the fast Fourier transform (FFT)

# Arguments:
- `audio_signal::Float`: the audio signal [number_samples, 1]
- `window_function::Float`: the window function [window_length, 1]
- `audio_mdct::Float`: the audio MDCT [number_frequencies, number_times]

# Example: Compute and display the MDCT as used in the AC-3 audio coding format
```
# Audio signal averaged over its channels and sample rate in Hz
Pkg.add("WAV")
using WAV
audio_signal, sample_rate = wavread("audio_file.wav");
audio_signal = mean(audio_signal, 2);

# Kaiser-Bessel-derived (KBD) window as used in the AC-3 audio coding format
window_length = 512;
alpha_value = 5;
include("z.jl")
window_function = z.kaiser(convert(Int64, window_length/2)+1, alpha_value*pi);
window_function2 = cumsum(window_function[1:convert(Int64, window_length/2)]);
window_function = sqrt.([window_function2; window_function2[convert(Int64, window_length/2):-1:1]]./sum(window_function));

# MDCT
audio_mdct = z.mdct(audio_signal, window_function);

# MDCT displayed in dB, s, and kHz
Pkg.add("Plots")
using Plots
plotly()
x_labels = [string(round(i*convert(Int64, window_length/2)/sample_rate, 2)) for i = 1:size(audio_mdct, 2)];
y_labels = [string(round(i*sample_rate/window_length/1000, 2)) for i = 1:size(audio_mdct, 1)];
heatmap(x_labels, y_labels, 20*log10.(abs.(audio_mdct)))
```
"""
function mdct(audio_signal, window_function)

    # Number of samples and window length
    number_samples = length(audio_signal);
    window_length = length(window_function);

    # Number of time frames
    number_times = ceil(Int64, 2*number_samples/window_length)+1;

    # Pre and post zero-padding of the signal
    audio_signal = [zeros(convert(Int64, window_length/2), 1); audio_signal;
    zeros(convert(Int64, (number_times+1)*window_length/2-number_samples), 1)];

    # Initialize the MDCT
    audio_mdct = zeros(convert(Int64, window_length/2), number_times);

    # Pre and post-processing arrays
    preprocessing_array = exp.(-im*pi/window_length*(0:window_length-1));
    postprocessing_array = exp.(-im*pi/window_length*(window_length/2+1)*(0.5:window_length/2-0.5));

    # Loop over the time frames
    for time_index = 1:number_times

        # Window the signal
        sample_index = convert(Int64, window_length/2)*(time_index-1);
        audio_segment = audio_signal[1+sample_index:window_length+sample_index].*window_function;

        # FFT of the audio segment after pre-processing
        audio_segment = fft(audio_segment.*preprocessing_array, 1);

        # Truncate to the first half before post-processing
        audio_mdct[:, time_index] = real(audio_segment[1:convert(Int64, window_length/2)].*postprocessing_array);

    end

    return audio_mdct

end

"""
audio_signal = z.imdct(audio_mdct, window_function);

    Compute the inverse modified discrete cosine transform (MDCT) using the fast Fourier transform (FFT)

# Arguments:
- `audio_mdct::Float`: the audio MDCT [number_frequencies, number_times]
- `window_function::Float`: the window function [window_length, 1]
- `audio_signal::Float`: the audio signal [number_samples, 1]

# Example: Verify that the MDCT is perfectly invertible
```
# Import modules
Pkg.add("WAV")
using WAV
audio_signal, sample_rate = wavread("audio_file.wav");
audio_signal = mean(audio_signal, 2);

# MDCT with a slope function as used in the Vorbis audio coding format
window_length = 2048;
window_function = sin.(pi/2*(sin.(pi/window_length*(0.5:window_length-0.5)).^2));
include("z.jl")
audio_mdct = z.mdct(audio_signal, window_function);

# Inverse MDCT and error signal
audio_signal2 = z.imdct(audio_mdct, window_function);
audio_signal2 = audio_signal2[1:length(audio_signal)];
error_signal = audio_signal-audio_signal2;

# Original, resynthesized, and error signals displayed in s
Pkg.add("Plots")
using Plots
plotly()
time_signal = (1:size(audio_signal, 1))/sample_rate;
audio_plot = plot(time_signal, audio_signal, xlabel="Time (s)", title="Original Signal");
audio2_plot = plot(time_signal, audio_signal2, xlabel="Time (s)", title="Resynthesized Signal");
error_plot = plot(time_signal, error_signal, xlabel="Time (s)", title="Error Signal");
plot(audio_plot, audio2_plot, error_plot, layout=(3,1), legend=false)
```
"""
function imdct(audio_mdct, window_function)

    # Number of frequency channels and time frames
    number_frequencies, number_times = size(audio_mdct);

    # Number of samples for the signal
    number_samples = number_frequencies*(number_times+1);

    # Initialize the audio signal
    audio_signal = zeros(number_samples, 1);

    # Pre and post-processing arrays
    preprocessing_array = exp.(-im*pi/(2*number_frequencies)*(number_frequencies+1)*(0:number_frequencies-1));
    postprocessing_array = exp.(-im*pi/(2*number_frequencies)*
    (0.5+number_frequencies/2:2*number_frequencies+number_frequencies/2-0.5))/number_frequencies;

    # FFT of the frames after pre-processing

    audio_mdct = fft([audio_mdct.*preprocessing_array;zeros(number_frequencies,number_times)], 1);

    # Apply the window to the frames after post-processing
    audio_mdct = 2*real(audio_mdct.*postprocessing_array).*window_function;

    # Loop over the time frames
    for time_index = 1:number_times

        # Recover the signal thanks to the time-domain aliasing cancellation (TDAC) principle
        sample_index = (time_index-1)*number_frequencies+1;
        audio_signal[sample_index:sample_index+2*number_frequencies-1] =
        audio_signal[sample_index:sample_index+2*number_frequencies-1]+audio_mdct[:, time_index];

    end

    # Remove the pre and post zero-padding
    audio_signal = audio_signal[number_frequencies+1:end-number_frequencies];

    return audio_signal

end

"Compute the Hamming window"
function hamming(window_length, window_sampling="symmetric")

    if window_sampling == "symmetric"
        window_function = 0.54 - 0.46*cos.(2*pi*(0:window_length-1)/(window_length-1));
    elseif window_sampling == "periodic"
        window_function = 0.54 - 0.46*cos.(2*pi*(0:window_length-1)/window_length);
    end

end

"Compute the Kaiser window"
function kaiser(window_length, alpha_value)

    window_function = zeros(window_length, 1);
    for window_index = 1:window_length
        window_function[window_index] = besseli(0, alpha_value*sqrt(1-(2*(window_index-1)/(window_length-1)-1).^2));
    end
    window_function = window_function/besseli(0, alpha_value);

end

end
