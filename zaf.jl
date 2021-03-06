"""
This Julia module implements a number of functions for audio signal analysis.

# Functions:
    stft - Compute the short-time Fourier transform (STFT).
    istft - Compute the inverse STFT.
    melfilterbank - Compute the mel filterbank.
    melspectrogram - Compute the mel spectrogram using a mel filterbank.
    mfcc - Compute the mel-frequency cepstral coefficients (MFCCs) using a mel filterbank.
    cqtkernel - Compute the constant-Q transform (CQT) kernel.
    cqtspectrogram - Compute the CQT spectrogram using a CQT kernel.
    cqtchromagram - Compute the CQT chromagram using a CQT kernel.
    dct - Compute the discrete cosine transform (DCT) using the fast Fourier transform (FFT).
    dst - Compute the discrete sine transform (DST) using the FFT.
    mdct - Compute the modified discrete cosine transform (MDCT) using the FFT.
    imdct - Compute the inverse MDCT using the FFT.

# Other:
    hamming - Compute the Hamming window.
    sigplot - Plot a signal in seconds.
    specshow - Display an spectrogram in dB, seconds, and Hz.
    melspecshow - Display a mel spectrogram in dB, seconds, and Hz.
    mfccshow - Display MFCCs in seconds.
    cqtspecshow - Display a CQT spectrogram in dB, seconds, and Hz.
    cqtchromshow - Display a CQT chromagram in seconds.

# Author:
    Zafar Rafii
    zafarrafii@gmail.com
    http://zafarrafii.com
    https://github.com/zafarrafii
    https://www.linkedin.com/in/zafarrafii/
    05/11/20
"""
module zaf

using FFTW, SparseArrays, Plots

export stft, istft, melfilterbank, melspectrogram, mfcc, cqtkernel, 
cqtspectrogram, cqtchromagram, dct, dst, mdct, imdct


"""
    audio_stft = stft(audio_signal, window_function, step_length)

Compute the short-time Fourier transform (STFT).

# Arguments:
- `audio_signal::Float`: the audio signal (number_samples,).
- `window_function::Float`: the window function (window_length,).
- `step_length::Integer`: the step length in samples.
- `audio_stft::Complex`: the audio STFT (window_length, number_frames).

# Example: Compute the spectrogram from an audio file.
```
# Load the needed modules
include("./zaf.jl")
using .zaf
using WAV
using Statistics
using Plots

# Read the audio signal with its sampling frequency in Hz, and average it over its channels
audio_signal, sampling_frequency = wavread("audio_file.wav")
audio_signal = mean(audio_signal, dims=2)

# Set the window duration in seconds (audio is stationary around 40 milliseconds)
window_duration = 0.04;

# Derive the window length in samples (use powers of 2 for faster FFT and constant overlap-add (COLA))
window_length = nextpow(2, ceil(Int, window_duration*sampling_frequency))

# Compute the window function (periodic Hamming window for COLA)
window_function = zaf.hamming(window_length, "periodic")

# Set the step length in samples (half of the window length for COLA)
step_length = convert(Int, window_length/2)

# Compute the STFT
audio_stft = zaf.stft(audio_signal, window_function, step_length)

# Derive the magnitude spectrogram (without the DC component and the mirrored frequencies)
audio_spectrogram = abs.(audio_stft[2:convert(Int, window_length/2)+1, :])

# Magnitude spectrogram (without the DC component and the mirrored frequencies)
audio_spectrogram = abs.(audio_stft[2:convert(Int, window_length/2)+1, :])

# Display the spectrogram in dB, seconds, and Hz
number_samples = length(audio_signal)
xtick_step = 1
ytick_step = 1000
plot_object = zaf.specshow(audio_spectrogram, number_samples, sampling_frequency, xtick_step, ytick_step)
heatmap!(title = "Spectrogram (dB)", size = (990, 600))
```
"""
function stft(audio_signal, window_function, step_length)

    # Get the number of samples and the window length in samples
    number_samples = length(audio_signal)
    window_length = length(window_function)

    # Derive the zero-padding length at the start and at the end of the signal to center the windows
    padding_length = floor(Int, window_length / 2)

    # Compute the number of time frames given the zero-padding at the start and at the end of the signal
    number_times =
        ceil(
            Int,
            ((number_samples + 2 * padding_length) - window_length) /
            step_length,
        ) + 1

    # Zero-pad the start and the end of the signal to center the windows
    audio_signal = [
        zeros(padding_length)
        audio_signal
        zeros(
            (
                number_times * step_length + (window_length - step_length) -
                padding_length
            ) - number_samples,
        )
    ]

    # Initialize the STFT
    audio_stft = zeros(window_length, number_times)

    # Loop over the time frames
    i = 0
    for j = 1:number_times

        # Window the signal
        audio_stft[:, j] = audio_signal[i+1:i+window_length] .* window_function
        i = i + step_length

    end

    # Compute the Fourier transform of the frames using the FFT
    audio_stft = fft(audio_stft, 1)

end

"""
    audio_istft = istft(audio_signal, window_function, step_length)

Compute the inverse short-time Fourier transform (STFT).

# Arguments:
- `audio_stft::Complex`: the audio STFT (window_length, number_frames).
- `window_function::Float`: the window function (window_length,).
- `step_length::Integer`: the step length in samples.
- `audio_signal::Float`: the audio signal (number_samples,).

# Example: Estimate the center and sides signals of a stereo audio file.
```
# Load the needed modules
include("./zaf.jl")
using .zaf
using WAV
using Plots

# Read the (stereo) audio signal with its sampling frequency in Hz
audio_signal, sampling_frequency = wavread("audio_file.wav")

# Set the parameters for the STFT
window_length = nextpow(2, ceil(Int, 0.04*sampling_frequency))
window_function = zaf.hamming(window_length, "periodic")
step_length = convert(Int, window_length/2)

# Compute the STFTs for the left and right channels
audio_stft1 = zaf.stft(audio_signal[:,1], window_function, step_length)
audio_stft2 = zaf.stft(audio_signal[:,2], window_function, step_length)

# Derive the magnitude spectrograms (with DC component) for the left and right channels
audio_spectrogram1 = abs.(audio_stft1[1:convert(Int, window_length/2)+1, :])
audio_spectrogram2 = abs.(audio_stft2[1:convert(Int, window_length/2)+1, :])

# Estimate the time-frequency masks for the left and right channels for the center
center_mask1 = min.(audio_spectrogram1, audio_spectrogram2)./audio_spectrogram1
center_mask2 = min.(audio_spectrogram1, audio_spectrogram2)./audio_spectrogram2

# Derive the STFTs for the left and right channels for the center (with mirrored frequencies)
center_stft1 = [center_mask1; center_mask1[convert(Int, window_length/2):-1:2,:]] .* audio_stft1
center_stft2 = [center_mask2; center_mask2[convert(Int, window_length/2):-1:2,:]] .* audio_stft2

# Synthesize the signals for the left and right channels for the center
center_signal1 = zaf.istft(center_stft1, window_function, step_length)
center_signal2 = zaf.istft(center_stft2, window_function, step_length)

# Derive the final stereo center and sides signals
center_signal = [center_signal1 center_signal2];
center_signal = center_signal[1:size(audio_signal, 1), :]
sides_signal = audio_signal-center_signal;

# Write the center and sides signals
wavwrite(center_signal, "center_signal.wav", Fs=sampling_frequency)
wavwrite(sides_signal, "sides_signal.wav", Fs=sampling_frequency)

# Display the original, center, and sides signals in seconds
xtick_step = 1
plot_object1 = zaf.sigplot(audio_signal, sampling_frequency, xtick_step)
plot!(ylims = (-1, 1), title = "Original signal")
plot_object2 = zaf.sigplot(center_signal, sampling_frequency, xtick_step)
plot!(ylims = (-1, 1), title = "Center signal")
plot_object3 = zaf.sigplot(sides_signal, sampling_frequency, xtick_step)
plot!(ylims = (-1, 1), title = "Sides signal")
plot(plot_object1, plot_object2, plot_object3, layout = (3, 1), size = (990, 600))
```
"""
function istft(audio_stft, window_function, step_length)

    # Get the window length in samples and the number of time frames
    window_length, number_times = size(audio_stft)

    # Compute the number of samples for the signal
    number_samples = number_times * step_length + (window_length - step_length)

    # Initialize the signal
    audio_signal = zeros(number_samples)

    # Compute the inverse Fourier transform of the frames and real part to ensure real values
    audio_stft = real(ifft(audio_stft, 1))

    # Loop over the time frames
    i = 0
    for j = 1:number_times

        # Perform a constant overlap-add (COLA) of the signal
        # (with proper window function and step length)
        audio_signal[i+1:i+window_length] =
            audio_signal[i+1:i+window_length] + audio_stft[:, j]
        i = i + step_length

    end

    # Remove the zero-padding at the start and at the end of the signal
    audio_signal =
        audio_signal[window_length-step_length+1:number_samples-(window_length-step_length)]

    # Normalize the signal by the gain introduced by the COLA (if any)
    audio_signal =
        audio_signal / sum(window_function[1:step_length:window_length])

end

"""
mel_filterbank = melfilterbank(sampling_frequency, window_length, number_filters)

    Compute the mel filterbank.

# Arguments:
- `sampling_frequency::Float`: the sampling frequency in Hz.
- `window_length::Integer`: the window length for the Fourier analysis in samples.
- `number_filters::Integer`: the number of mel filters.
- `mel_filterbank::Float`: the mel filterbank (sparse) (number_mels, number_frequencies).

# Example: Compute and display the mel filterbank.
```
# Load the needed modules
include("./zaf.jl")
using .zaf
using Plots

# Compute the mel filterbank using some parameters
sampling_frequency = 44100
window_length = nextpow(2, ceil(Int, 0.04*sampling_frequency))
number_mels = 128
mel_filterbank = zaf.melfilterbank(sampling_frequency, window_length, number_mels)

# Display the mel filterbank
heatmap(Array(mel_filterbank), fillcolor = :jet, legend = false, fmt = :png, size = (990, 300), 
    title = "Mel filterbank", xlabel = "Frequency index", ylabel = "Mel index")
```
"""
function melfilterbank(
    sampling_frequency,
    window_length,
    number_filters,
)

    # Compute the minimum and maximum mels
    minimum_mel = 2595 * log10(1 + (sampling_frequency / window_length) / 700)
    maximum_mel = 2595 * log10(1 + (sampling_frequency / 2) / 700)
    
    # Derive the width of the half-overlapping filters in the mel scale (constant)
    filter_width =
        2 * (maximum_mel - minimum_mel) / (number_filters + 1)

    # Compute the start and end indices of the filters in the mel scale (linearly spaced)
    filter_indices = [minimum_mel:filter_width/2:maximum_mel;]

    # Derive the indices of the filters in the linear frequency scale (log spaced)
    filter_indices =
        round.(
            Int,
            700 * (10 .^ (filter_indices / 2595) .- 1) * window_length /
            sampling_frequency,
        )

    # Initialize the mel filterbank
    mel_filterbank = zeros(number_filters, convert(Int, window_length / 2))

    # Loop over the filters
    for i = 1:number_filters

        # Compute the left and right sides of the triangular filters
        # (this is more accurate than creating triangular filters directly)
        mel_filterbank[i, filter_indices[i]:filter_indices[i+1]] = range(
            0,
            stop = 1,
            length = filter_indices[i+1] - filter_indices[i] + 1,
        )
        mel_filterbank[i, filter_indices[i+1]:filter_indices[i+2]] = range(
            1,
            stop = 0,
            length = filter_indices[i+2] - filter_indices[i+1] + 1,
        )

    end

    # Make the mel filterbank sparse
    mel_filterbank = sparse(mel_filterbank)

end

"""
melspectrogram(audio_signal, window_function, step_length, mel_filterbank)

    Compute the mel spectrogram using a mel filterbank.

# Arguments:
- `audio_signal::Float`: the audio signal (number_samples,).
- `window_function::Float`: the window function (window_length,).
- `step_length::Integer`: the step length in samples.
- `mel_filterbank::Float`: the mel filterbank (number_mels, number_frequencies).
- `mel_spectrogram::Float`: the mel spectrogram (number_mels, number_times).

# Example: Compute and display the mel spectrogram.
```
# Load the needed modules
include("./zaf.jl")
using .zaf
using WAV
using Statistics
using Plots

# Read the audio signal with its sampling frequency in Hz, and average it over its channels
audio_signal, sampling_frequency = wavread("audio_file.wav")
audio_signal = mean(audio_signal, dims=2)

# Set the parameters for the Fourier analysis
window_length = nextpow(2, ceil(Int, 0.04*sampling_frequency))
window_function = zaf.hamming(window_length, "periodic")
step_length = convert(Int, window_length/2)

# Compute the mel filterbank
number_mels = 128
mel_filterbank = zaf.melfilterbank(sampling_frequency, window_length, number_mels)

# Compute the mel spectrogram using the filterbank
mel_spectrogram = zaf.melspectrogram(audio_signal, window_function, step_length, mel_filterbank)

# Display the mel spectrogram in dB, seconds, and Hz
number_samples = length(audio_signal)
xtick_step = 1
plot_object = zaf.melspecshow(mel_spectrogram, number_samples, sampling_frequency, window_length, xtick_step)
heatmap!(title = "Mel spectrogram (dB)", size = (990, 600))
```
"""
function melspectrogram(
    audio_signal,
    window_function,
    step_length,
    mel_filterbank,
)

    # Compute the magnitude spectrogram (without the DC component and the mirrored frequencies)
    audio_stft = zaf.stft(audio_signal, window_function, step_length)
    audio_spectrogram = abs.(audio_stft[2:convert(Int, length(window_function)/2)+1, :])
    
    # Compute the mel spectrogram by using the filterbank
    mel_spectrogram = mel_filterbank * audio_spectrogram

end

"""
audio_mfcc = mfcc(audio_signal, window_function, step_length, mel_filterbank, number_coefficients)

    Compute the mel-frequency cepstral coefficients (MFCCs).

# Arguments:
- `audio_signal::Float`: the audio signal (number_samples,).
- `window_function::Float`: the window function (window_length,).
- `step_length::Integer`: the step length in samples.
- `mel_filterbank::Float`: the mel filterbank (number_mels, number_frequencies).
- `number_coefficients::Integer`: the number of coefficients (without the 0th coefficient).
- `audio_mfcc::Float`: the audio MFCCs (number_coefficients, number_times).

# Example: Compute and display the MFCCs, delta MFCCs, and delta-detla MFCCs.
```
# Load the needed modules
include("./zaf.jl")
using .zaf
using WAV
using Statistics
using Plots

# Read the audio signal with its sampling frequency in Hz, and average it over its channels
audio_signal, sampling_frequency = wavread("audio_file.wav")
audio_signal = mean(audio_signal, dims=2)

# Set the parameters for the Fourier analysis
window_length = nextpow(2, ceil(Int, 0.04*sampling_frequency))
window_function = zaf.hamming(window_length, "periodic")
step_length = convert(Int, window_length/2)

# Compute the mel filterbank
number_mels = 40
mel_filterbank = zaf.melfilterbank(sampling_frequency, window_length, number_mels)

# Compute the MFCCs using the filterbank
number_coefficients = 20
audio_mfcc = zaf.mfcc(audio_signal, window_function, step_length, mel_filterbank, number_coefficients)

# Compute the delta and delta-delta MFCCs
audio_dmfcc = diff(audio_mfcc, dims=2)
audio_ddmfcc = diff(audio_dmfcc, dims=2)

# Display the MFCCs, delta MFCCs, and delta-delta MFCCs in seconds
xtick_step = 1
number_samples = length(audio_signal)
plot_object1 = zaf.mfccshow(audio_mfcc, number_samples, sampling_frequency, xtick_step); plot!(title = "MFCCs")
plot_object2 = zaf.mfccshow(audio_dmfcc, number_samples, sampling_frequency, xtick_step); plot!(title = "Delta MFCCs")
plot_object3 = zaf.mfccshow(audio_ddmfcc, number_samples, sampling_frequency, xtick_step); plot!(title = "Delta MFCCs")
plot(plot_object1, plot_object2, plot_object3, layout = (3, 1), size = (990, 600))
```
"""
function mfcc(
    audio_signal,
    window_function,
    step_length,
    mel_filterbank, 
    number_coefficients,
)

    # Compute the power spectrogram (without the DC component and the mirrored frequencies)
    audio_stft = zaf.stft(audio_signal, window_function, step_length)
    audio_spectrogram = 2 .^ abs.(audio_stft[2:convert(Int, length(window_function)/2)+1, :])

    # Compute the discrete cosine transform of the log magnitude spectrogram
    # mapped onto the mel scale using the filter bank
    audio_mfcc = FFTW.dct(log.(mel_filterbank * audio_spectrogram .+ eps()), 1)

    # Keep only the first coefficients (without the 0th)
    audio_mfcc = audio_mfcc[2:number_coefficients+1, :]

end

"""
    cqt_kernel = cqtkernel(sampling_frequency, frequency_resolution, minimum_frequency, maximum_frequency)

Compute the constant-Q transform (CQT) kernel.

# Arguments:
- `sampling_frequency::Float` the sampling frequency in Hz.
- `frequency_resolution::Integer` the frequency resolution in number of frequency channels per semitone.
- `minimum_frequency::Float`: the minimum frequency in Hz.
- `maximum_frequency::Float`: the maximum frequency in Hz.
- `cqt_kernel::Complex`: the CQT kernel (sparse) (number_frequencies, fft_length).

# Example: Compute and display the CQT kernel.
```
# Load the needed modules
include("./zaf.jl")
using .zaf
using Plots

# Set the parameters for the CQT kernel
sampling_frequency = 44100
octave_resolution = 24
minimum_frequency = 55
maximum_frequency = sampling_frequency/2

# Compute the CQT kernel
cqt_kernel = zaf.cqtkernel(sampling_frequency, octave_resolution, minimum_frequency, maximum_frequency)

# Display the magnitude CQT kernel
heatmap(abs.(Array(cqt_kernel)), fillcolor = :jet, legend = false, fmt = :png, size = (990, 300),
    title = "Magnitude CQT kernel", xlabel = "FFT index", ylabel = "CQT index")
```
"""
function cqtkernel(
    sampling_frequency,
    octave_resolution,
    minimum_frequency,
    maximum_frequency,
)

    # Compute the constant ratio of frequency to resolution (= fk/(fk+1-fk))
    quality_factor = 1 / (2^(1 / octave_resolution) - 1)

    # Compute the number of frequency channels for the CQT
    number_frequencies = round(
        Int,
        octave_resolution * log2(maximum_frequency / minimum_frequency),
    )

    # Compute the window length for the FFT (= longest window for the minimum frequency)
    fft_length = nextpow(
        2,
        ceil(Int, quality_factor * sampling_frequency / minimum_frequency),
    )

    # Initialize the (complex) CQT kernel
    cqt_kernel = zeros(ComplexF64, number_frequencies, fft_length)

    # Loop over the frequency channels
    for i = 1:number_frequencies

        # Derive the frequency value in Hz
        frequency_value = minimum_frequency * 2^((i - 1) / octave_resolution)

        # Compute the window length in samples
        # (nearest odd value to center the temporal kernel on 0)
        window_length =
            2 * round(
                Int,
                quality_factor * sampling_frequency / frequency_value / 2,
            ) + 1

        # Compute the temporal kernel for the current frequency (odd and symmetric)
        temporal_kernel =
            hamming(window_length, "symmetric") .*
            exp.(
                2 *
                pi *
                im *
                quality_factor *
                (-(window_length - 1)/2:(window_length-1)/2) / window_length,
            ) / window_length

        # Derive the pad width to center the temporal kernels
        pad_width = convert(Int, (fft_length - window_length + 1) / 2)

        # Save the current temporal kernel at the center
        # (the zero-padded temporal kernels are not perfectly symmetric anymore because of the even length here)
        cqt_kernel[i, pad_width+1:pad_width+window_length] = temporal_kernel

    end

    # Derive the spectral kernels by taking the FFT of the temporal kernels
    # (the spectral kernels are almost real because the temporal kernels are almost symmetric)
    cqt_kernel = fft(cqt_kernel, 2)

    # Make the CQT kernel sparser by zeroing magnitudes below a threshold
    cqt_kernel[abs.(cqt_kernel).<0.01] .= 0

    # Make the CQT kernel sparse
    cqt_kernel = sparse(cqt_kernel)

    # Get the final CQT kernel by using Parseval's theorem
    cqt_kernel = conj.(cqt_kernel) / fft_length

end

"""
    audio_spectrogram = cqtspectrogram(audio_signal, sampling_frequency, time_resolution, cqt_kernel)

Compute the constant-Q transform (CQT) spectrogram using a kernel.

# Arguments:
- `audio_signal::Float`: the audio signal (number_samples,).
- `sample_rate::Float`: the sample rate in Hz.
- `time_resolution::Float`: the number of time frames per second.
- `cqt_kernel::Complex`: the CQT kernel (number_frequencies, fft_length).
- `cqt_spectrogram::Float`: the CQT spectrogram (number_frequencies, number_times).

# Example: Compute and display the CQT spectrogram.
```
# Load the needed modules
include("./zaf.jl")
using .zaf
using WAV
using Statistics
using Plots

# Read the audio signal with its sampling frequency in Hz, and average it over its channels
audio_signal, sampling_frequency = wavread("audio_file.wav")
audio_signal = mean(audio_signal, dims=2)

# Compute the CQT kernel
octave_resolution = 24
minimum_frequency = 55
maximum_frequency = 3520
cqt_kernel = zaf.cqtkernel(sampling_frequency, octave_resolution, minimum_frequency, maximum_frequency)

# Compute the CQT spectrogram using the kernel
time_resolution = 25
audio_spectrogram = zaf.cqtspectrogram(audio_signal, sampling_frequency, time_resolution, cqt_kernel)

# Display the CQT spectrogram in dB, seconds, and Hz
xtick_step = 1
plot_object = zaf.cqtspecshow(cqt_spectrogram, time_resolution, octave_resolution, minimum_frequency, xtick_step)
heatmap!(title = "CQT spectrogram (dB)", size = (990, 600))
```
"""
function cqtspectrogram(
    audio_signal,
    sampling_frequency,
    time_resolution,
    cqt_kernel,
)

    # Derive the number of time samples per time frame
    step_length = round(Int, sampling_frequency / time_resolution)

    # Compute the number of time frames
    number_times = floor(Int, length(audio_signal) / step_length)

    # Get the number of frequency channels and the FFT length
    number_frequencies, fft_length = size(cqt_kernel)

    # Zero-pad the signal to center the CQT
    audio_signal = [
        zeros(ceil(Int, (fft_length - step_length) / 2))
        audio_signal
        zeros(floor(Int, (fft_length - step_length) / 2))
    ]

    # Initialize the spectrogram
    audio_spectrogram = zeros(number_frequencies, number_times)

    # Loop over the time frames
    i = 0
    for j = 1:number_times

        # Compute the magnitude CQT using the kernel
        audio_spectrogram[:, j] =
            abs.(cqt_kernel * fft(audio_signal[i+1:i+fft_length]))
        i = i + step_length

    end

    # Return the output explicitly
    return audio_spectrogram

end

"""
    cqt_chromagram = cqtchromagram(audio_signal, sampling_frequency, time_resolution, octave_resolution, cqt_kernel)

Compute the constant-Q transform (CQT) chromagram using a kernel.

# Arguments:
- `audio_signal::Float`: the audio signal (number_samples,).
- `sampling_frequency::Float`: the sample rate in Hz.
- `time_resolution::Float`: the number of time frames per second.
- `octave_resolution::Integer`: the number of frequency channels per octave.
- `cqt_kernel::Complex`: the CQT kernel (number_frequencies, fft_length).
- `cqt_chromagram::Float`: the CQT chromagram (number_chromas, number_times).

# Example: Compute and display the CQT chromagram.
```
# Load the needed modules
include("./zaf.jl")
using .zaf
using WAV
using Statistics
using Plots

# Read the audio signal with its sampling frequency in Hz, and average it over its channels
audio_signal, sampling_frequency = wavread("audio_file.wav")
audio_signal = mean(audio_signal, dims=2)

# Compute the CQT kernel
octave_resolution = 24
minimum_frequency = 55
maximum_frequency = 3520
cqt_kernel = zaf.cqtkernel(sampling_frequency, octave_resolution, minimum_frequency, maximum_frequency)

# Compute the CQT chromagram using the kernel
time_resolution = 25
cqt_chromagram = zaf.cqtchromagram(audio_signal, sampling_frequency, time_resolution, octave_resolution, cqt_kernel)

# Display the CQT chromagram in seconds
xtick_step = 1
plot_object = zaf.cqtchromshow(cqt_chromagram, time_resolution, xtick_step)
heatmap!(title = "CQT chromagram", size = (990, 300))
```
"""
function cqtchromagram(
    audio_signal,
    sampling_frequency,
    time_resolution,
    octave_resolution,
    cqt_kernel,
)

    # Compute the CQT spectrogram
    cqt_spectrogram = cqtspectrogram(
        audio_signal,
        sampling_frequency,
        time_resolution,
        cqt_kernel,
    )

    # Get the number of frequency channels and time frames
    number_frequencies, number_times = size(cqt_spectrogram)

    # Initialize the chromagram
    cqt_chromagram = zeros(octave_resolution, number_times)

    # Loop over the chroma channels
    for i = 1:octave_resolution

        # Sum the energy of the frequency channels for every chroma
        cqt_chromagram[i, :] = sum(
            cqt_spectrogram[i:octave_resolution:number_frequencies, :],
            dims = 1,
        )

    end

    # Return the output explicitly
    return cqt_chromagram

end

"""
audio_dct = dct(audio_signal, dct_type)

    Compute the discrete cosine transform (DCT) using the fast Fourier transform (FFT).

# Arguments:
- `audio_signal::Float`: the audio signal (window_length,).
- dct_type::Integer`: the DCT type (1, 2, 3, or 4).
- audio_dct::Float`: the audio DCT (number_frequencies,).

# Example: Compute the 4 different DCTs and compare them to FFTW's DCTs.
```
# Load the needed modules
include("./zaf.jl")
using .zaf
using WAV
using Statistics
using FFTW
using Plots

# Read the audio signal with its sampling frequency in Hz, and average it over its channels
audio_signal, sampling_frequency = wavread("audio_file.wav")
audio_signal = mean(audio_signal, dims=2)

# Get an audio segment for a given window length
window_length = 1024
audio_segment = audio_signal[1:window_length]

# Compute the DCT-I, II, III, and IV
audio_dct1 = zaf.dct(audio_segment, 1)
audio_dct2 = zaf.dct(audio_segment, 2)
audio_dct3 = zaf.dct(audio_segment, 3)
audio_dct4 = zaf.dct(audio_segment, 4)

# Compute FFTW's DCT-I, II, III, and IV (orthogonalized)
audio_segment2 = copy(audio_segment)
audio_segment2[[1, end]] = audio_segment2[[1, end]]*sqrt(2)
fftw_dct1 = FFTW.r2r(audio_segment2, FFTW.REDFT00)
fftw_dct1[[1, window_length]] = fftw_dct1[[1, window_length]]/sqrt(2)
fftw_dct1 = fftw_dct1/2 * sqrt(2/(window_length - 1))

fftw_dct2 = FFTW.r2r(audio_segment, FFTW.REDFT10)
fftw_dct2[1] = fftw_dct2[1]/sqrt(2)
fftw_dct2 = fftw_dct2/2 * sqrt(2/window_length)

audio_segment2 = copy(audio_segment)
audio_segment2[1] = audio_segment2[1]*sqrt(2)
fftw_dct3 = FFTW.r2r(audio_segment2, FFTW.REDFT01)
fftw_dct3 = fftw_dct3/2 * sqrt(2/window_length)

fftw_dct4 = FFTW.r2r(audio_segment, FFTW.REDFT11)
fftw_dct4 = fftw_dct4/2 * sqrt(2/window_length)

# Plot the DCT-I, II, III, and IV, FFTW's versions, and their differences (using yformatter because of precision issue in ticks)
dct1_plot = plot(audio_dct1, title = "DCT-I")
dct2_plot = plot(audio_dct2, title = "DCT-II")
dct3_plot = plot(audio_dct3, title = "DCT-III")
dct4_plot = plot(audio_dct4, title = "DCT-IV")
dct1_plot2 = plot(fftw_dct1, title = "FFTW's DCT-I")
dct2_plot2 = plot(fftw_dct2, title = "FFTW's DCT-II")
dct3_plot2 = plot(fftw_dct3, title = "FFTW's DCT-III")
dct4_plot2 = plot(fftw_dct4, title = "FFTW's DCT-IV")
diff1_plot = plot(audio_dct1-fftw_dct1, title = "DCT-I - FFTW's DCT-I")
diff2_plot = plot(audio_dct2-fftw_dct2, title = "DCT-II - FFTW's DCT-II", yformatter = y->string(convert(Int, round(y/1e-16)),"x10⁻¹⁶"))
diff3_plot = plot(audio_dct3-fftw_dct3, title = "DCT-III - FFTW's DCT-III", yformatter = y->string(convert(Int, round(y/1e-16)),"x10⁻¹⁶"))
diff4_plot = plot(audio_dct4-fftw_dct4, title = "DCT-IV - FFTW's DCT-IV", yformatter = y->string(convert(Int, round(y/1e-16)),"x10⁻¹⁶"))
plot(dct1_plot, dct2_plot, dct3_plot, dct4_plot, dct1_plot2, dct2_plot2, dct3_plot2, dct4_plot2,
    diff1_plot, diff2_plot, diff3_plot, diff4_plot, xlims = (0, window_length), legend = false, titlefont = 10,
    layout = (3,4), size = (990, 600), fmt = :png)
```
"""
function dct(audio_signal, dct_type)

    # Check if the DCT type is I, II, III, or IV
    if dct_type == 1

        # Get the number of samples
        window_length = length(audio_signal)

        # Pre-process the signal to make the DCT-I matrix orthogonal
        # (copy the signal to avoid modifying it outside of the function)
        audio_signal = copy(audio_signal)
        audio_signal[[1, end]] = audio_signal[[1, end]] * sqrt(2)

        # Compute the DCT-I using the FFT
        audio_dct = [audio_signal; audio_signal[end-1:-1:2]]
        audio_dct = fft(audio_dct, 1)
        audio_dct = real(audio_dct[1:window_length]) / 2

        # Post-process the results to make the DCT-I matrix orthogonal
        audio_dct[[1, end]] = audio_dct[[1, end]] / sqrt(2)
        audio_dct = audio_dct * sqrt(2 / (window_length - 1))

    elseif dct_type == 2

        # Get the number of samples
        window_length = length(audio_signal)

        # Compute the DCT-II using the FFT
        audio_dct = zeros(4 * window_length)
        audio_dct[2:2:2*window_length] = audio_signal
        audio_dct[2*window_length+2:2:4*window_length] = audio_signal[end:-1:1]
        audio_dct = fft(audio_dct, 1)
        audio_dct = real(audio_dct[1:window_length]) / 2

        # Post-process the results to make the DCT-II matrix orthogonal
        audio_dct[1] = audio_dct[1] / sqrt(2)
        audio_dct = audio_dct * sqrt(2 / window_length)

    elseif dct_type == 3

        # Get the number of samples
        window_length = length(audio_signal)

        # Pre-process the signal to make the DCT-III matrix orthogonal
        # (copy the signal to avoid modifying it outside of the function)
        audio_signal = copy(audio_signal)
        audio_signal[1] = audio_signal[1] * sqrt(2)

        # Compute the DCT-III using the FFT
        audio_dct = zeros(4 * window_length)
        audio_dct[1:window_length] = audio_signal
        audio_dct[window_length+2:2*window_length+1] = -audio_signal[end:-1:1]
        audio_dct[2*window_length+2:3*window_length] = -audio_signal[2:end]
        audio_dct[3*window_length+2:4*window_length] = audio_signal[end:-1:2]
        audio_dct = fft(audio_dct, 1)
        audio_dct = real(audio_dct[2:2:2*window_length]) / 4

        # Post-process the results to make the DCT-III matrix orthogonal
        audio_dct = audio_dct * sqrt(2 / window_length)

    elseif dct_type == 4

        # Get the number of samples
        window_length = length(audio_signal)

        # Compute the DCT-IV using the FFT
        audio_dct = zeros(8 * window_length)
        audio_dct[2:2:2*window_length] = audio_signal
        audio_dct[2*window_length+2:2:4*window_length] = -audio_signal[end:-1:1]
        audio_dct[4*window_length+2:2:6*window_length] = -audio_signal
        audio_dct[6*window_length+2:2:8*window_length] = audio_signal[end:-1:1]
        audio_dct = fft(audio_dct, 1)
        audio_dct = real(audio_dct[2:2:2*window_length]) / 4

        # Post-process the results to make the DCT-IV matrix orthogonal
        audio_dct = audio_dct * sqrt(2 / window_length)

    end

end

"""
audio_dst = dst(audio_signal, dst_type)

    Compute the discrete sine transform (DST) using the fast Fourier transform (FFT).

# Arguments:
- `audio_signal::Float`: the audio signal (window_length,).
- `dst_type::Integer`: the DST type (1, 2, 3, or 4).
- `audio_dst::Float`: the audio DST (number_frequencies,)

# Example: Compute the 4 different DSTs and compare them to their respective inverses.
```
# Load the needed modules
include("./zaf.jl")
using .zaf
using WAV
using Statistics
using Plots

# Read the audio signal with its sampling frequency in Hz, and average it over its channels
audio_signal, sampling_frequency = wavread("audio_file.wav")
audio_signal = mean(audio_signal, dims=2)

# Get an audio segment for a given window length
window_length = 1024
audio_segment = audio_signal[1:window_length]

# Compute the DST-I, II, III, and IV
audio_dst1 = zaf.dst(audio_segment, 1)
audio_dst2 = zaf.dst(audio_segment, 2)
audio_dst3 = zaf.dst(audio_segment, 3)
audio_dst4 = zaf.dst(audio_segment, 4)

# Compute their respective inverses, i.e., DST-I, II, III, and IV
audio_idst1 = zaf.dst(audio_dst1, 1)
audio_idst2 = zaf.dst(audio_dst2, 3)
audio_idst3 = zaf.dst(audio_dst3, 2)
audio_idst4 = zaf.dst(audio_dst4, 4)

# Plot the DST-I, II, III, and IV, their respective inverses, and their differences with the original audio segment
dst1_plot = plot(audio_dst1, title = "DST-I")
dst2_plot = plot(audio_dst2, title = "DST-II")
dst3_plot = plot(audio_dst3, title = "DST-III")
dst4_plot = plot(audio_dst4, title = "DST-IV")
idst1_plot2 = plot(audio_idst1, title = "Inverse DST-I (DST-I)")
idst2_plot2 = plot(audio_idst2, title = "Inverse DST-II (DST-III)")
idst3_plot2 = plot(audio_idst3, title = "Inverse DST-III (DST-II)")
idst4_plot2 = plot(audio_idst4, title = "Inverse DST-IV (DST-IV)")
diff1_plot = plot(audio_idst1-audio_segment, title = "Inverse DST-I - audio segment")
diff2_plot = plot(audio_idst2-audio_segment, title = "Inverse DST-II - audio segment")
diff3_plot = plot(audio_idst3-audio_segment, title = "Inverse DST-III - audio segment")
diff4_plot = plot(audio_idst4-audio_segment, title = "Inverse DST-IV - audio segment")
plot(dst1_plot, dst2_plot, dst3_plot, dst4_plot, idst1_plot2, idst2_plot2, idst3_plot2, idst4_plot2,
    diff1_plot, diff2_plot, diff3_plot, diff4_plot, xlims = (0, window_length), legend = false, titlefont = 10,
    layout = (3,4), size = (990, 600), fmt = :png)
```
"""
function dst(audio_signal, dst_type)

    # Check if the DST type is I, II, III, or IV
    if dst_type == 1

        # Get the number of samples
        window_length = length(audio_signal)

        # Compute the DST-I using the FFT
        audio_dst = zeros(2 * window_length + 2)
        audio_dst[2:window_length+1] = audio_signal
        audio_dst[window_length+3:end] = -audio_signal[end:-1:1]
        audio_dst = fft(audio_dst, 1)
        audio_dst = -imag(audio_dst[2:window_length+1]) / 2

        # Post-process the results to make the DST-I matrix orthogonal
        audio_dst = audio_dst * sqrt(2 / (window_length + 1))

    elseif dst_type == 2

        # Get the number of samples
        window_length = length(audio_signal)

        # Compute the DST-II using the FFT
        audio_dst = zeros(4 * window_length)
        audio_dst[2:2:2*window_length] = audio_signal
        audio_dst[2*window_length+2:2:4*window_length] = -audio_signal[end:-1:1]
        audio_dst = fft(audio_dst, 1)
        audio_dst = -imag(audio_dst[2:window_length+1]) / 2

        # Post-process the results to make the DST-II matrix orthogonal
        audio_dst[end] = audio_dst[end] / sqrt(2)
        audio_dst = audio_dst * sqrt(2 / window_length)

    elseif dst_type == 3

        # Get the number of samples
        window_length = length(audio_signal)

        # Pre-process the signal to make the DST-III matrix orthogonal
        # (copy the signal to avoid modifying it outside of the function)
        audio_signal = copy(audio_signal)
        audio_signal[end] = audio_signal[end] * sqrt(2)

        # Compute the DST-III using the FFT
        audio_dst = zeros(4 * window_length)
        audio_dst[2:window_length+1] = audio_signal
        audio_dst[window_length+2:2*window_length] = audio_signal[end-1:-1:1]
        audio_dst[2*window_length+2:3*window_length+1] = -audio_signal
        audio_dst[3*window_length+2:4*window_length] = -audio_signal[end-1:-1:1]
        audio_dst = fft(audio_dst, 1)
        audio_dst = -imag(audio_dst[2:2:2*window_length]) / 4

        # Post-processing to make the DST-III matrix orthogonal
        audio_dst = audio_dst * sqrt(2 / window_length)

    elseif dst_type == 4

        # Get the number of samples
        window_length = length(audio_signal)

        # Compute the DST-IV using the FFT
        audio_dst = zeros(8 * window_length)
        audio_dst[2:2:2*window_length] = audio_signal
        audio_dst[2*window_length+2:2:4*window_length] = audio_signal[end:-1:1]
        audio_dst[4*window_length+2:2:6*window_length] = -audio_signal
        audio_dst[6*window_length+2:2:8*window_length] = -audio_signal[end:-1:1]
        audio_dst = fft(audio_dst, 1)
        audio_dst = -imag(audio_dst[2:2:2*window_length]) / 4

        # Post-process the results to make the DST-IV matrix orthogonal
        audio_dst = audio_dst * sqrt(2 / window_length)

    end

end

"""
audio_mdct = mdct(audio_signal, window_function)

    Compute the modified discrete cosine transform (MDCT) using the fast Fourier transform (FFT).

# Arguments:
- `audio_signal::Float`: the audio signal (number_samples,).
- `window_function::Float`: the window function (window_length,).
- `audio_mdct::Float`: the audio MDCT (number_frequencies, number_times).

# Example: Compute and display the MDCT as used in the AC-3 audio coding format.
```
# Load the needed modules
include("./zaf.jl")
using .zaf
using WAV
using Statistics
using Plots

# Add and use the DSP package to get the kaiser window function
using Pkg
Pkg.add("DSP")
using DSP

# Read the audio signal with its sampling frequency in Hz, and average it over its channels
audio_signal, sampling_frequency = wavread("audio_file.wav")
audio_signal = mean(audio_signal, dims=2)

# Compute the Kaiser-Bessel-derived (KBD) window as used in the AC-3 audio coding format
window_length = 512
alpha_value = 5
window_function = kaiser(convert(Int, window_length/2)+1, alpha_value*pi)
window_function2 = cumsum(window_function[1:convert(Int, window_length/2)])
window_function = sqrt.([window_function2; window_function2[convert(Int, window_length/2):-1:1]]./sum(window_function))

# Compute the MDCT
audio_mdct = zaf.mdct(audio_signal, window_function)

# Display the MDCT in dB, seconds, and Hz
xtick_step = 1
ytick_step = 1000
plot_object = zaf.specshow(abs.(audio_mdct), length(audio_signal), sampling_frequency, xtick_step, ytick_step)
heatmap!(title = "MDCT (dB)", size = (990, 600))
```
"""
function mdct(audio_signal, window_function)

    # Get the number of samples and the window length in samples
    number_samples = length(audio_signal)
    window_length = length(window_function)

    # Derive the step length and the number of frequencies (for clarity)
    step_length = convert(Int, window_length / 2)
    number_frequencies = convert(Int, window_length / 2)

    # Derive the number of time frames
    number_times = ceil(Int, number_samples / step_length) + 1

    # Zero-pad the start and the end of the signal to center the windows
    audio_signal = [
        zeros(step_length)
        audio_signal
        zeros((number_times + 1) * step_length - number_samples)
    ]

    # Initialize the MDCT
    audio_mdct = zeros(number_frequencies, number_times)

    # Prepare the pre-processing and post-processing arrays
    preprocessing_array = exp.(-im * pi / window_length * (0:window_length-1))
    postprocessing_array =
        exp.(
            -im * pi / window_length *
            (window_length / 2 + 1) *
            (0.5:window_length/2-0.5),
        )

    # Loop over the time frames
    # (do the pre and post-processing, and take the FFT in the loop to avoid storing twice longer frames)
    i = 0
    for j = 1:number_times

        # Window the signal
        audio_segment = audio_signal[i+1:i+window_length] .* window_function
        i = i + step_length

        # Compute the Fourier transform of the windowed segment using the FFT after pre-processing
        audio_segment = fft(audio_segment .* preprocessing_array)

        # Truncate to the first half before post-processing (and take the real to ensure real values)
        audio_mdct[:, j] =
            real(audio_segment[1:number_frequencies] .* postprocessing_array)

    end

    # Return the output explicitly as it is not clear here what the last value is
    return audio_mdct

end

"""
audio_signal = imdct(audio_mdct, window_function)

    Compute the inverse modified discrete cosine transform (MDCT) using the fast Fourier transform (FFT).

# Arguments:
- `audio_mdct::Float`: the audio MDCT (number_frequencies, number_times).
- `window_function::Float`: the window function (window_length,).
- `audio_signal::Float`: the audio signal (number_samples,).

# Example: Verify that the MDCT is perfectly invertible.
```
# Load the needed modules
include("./zaf.jl")
using .zaf
using WAV
using Statistics
using Plots

# Read the audio signal with its sampling frequency in Hz, and average it over its channels
audio_signal, sampling_frequency = wavread("audio_file.wav")
audio_signal = mean(audio_signal, dims=2)

# Compute the MDCT with a slope function as used in the Vorbis audio coding format
window_length = 2048
window_function = sin.(pi/2*(sin.(pi/window_length*(0.5:window_length-0.5)).^2))
audio_mdct = zaf.mdct(audio_signal, window_function)

# Compute the inverse MDCT
audio_signal2 = zaf.imdct(audio_mdct, window_function)
audio_signal2 = audio_signal2[1:length(audio_signal)]

# Compute the differences between the original signal and the resynthesized one
audio_differences = audio_signal-audio_signal2
y_max = maximum(abs.(audio_differences))

# Display the original and resynthesized signals, and their differences in seconds
xtick_step = 1
plot_object1 = zaf.sigplot(audio_signal, sampling_frequency, xtick_step)
plot!(ylims = (-1, 1), title = "Original signal")
plot_object2 = zaf.sigplot(audio_signal2, sampling_frequency, xtick_step)
plot!(ylims = (-1, 1), title = "Resyntesized signal")
plot_object3 = zaf.sigplot(audio_differences, sampling_frequency, xtick_step)
plot!(ylims = (-y_max, y_max), title = "Original - resyntesized signal")
plot(plot_object1, plot_object2, plot_object3, layout = (3, 1), size = (990, 600))
```
"""
function imdct(audio_mdct, window_function)

    # Get the number of frequency channels and time frames
    number_frequencies, number_times = size(audio_mdct)

    # Derive the window length and the step length in samples (for clarity)
    window_length = 2 * number_frequencies
    step_length = number_frequencies

    # Derive the number of samples for the signal
    number_samples = step_length * (number_times + 1)

    # Initialize the audio signal
    audio_signal = zeros(number_samples)

    # Prepare the pre-processing and post-processing arrays
    preprocessing_array =
        exp.(
            -im * pi / (2 * number_frequencies) *
            (number_frequencies + 1) *
            (0:number_frequencies-1),
        )
    postprocessing_array =
        exp.(
            -im * pi / (2 * number_frequencies) * (
                0.5+number_frequencies/2:2*number_frequencies+number_frequencies/2-0.5
            ),
        ) / number_frequencies

    # Compute the Fourier transform of the frames using the FFT after pre-processing
    # (zero-pad to get twice the length)
    audio_mdct = fft(
        [
            audio_mdct .* preprocessing_array
            zeros(number_frequencies, number_times)
        ],
        1,
    )

    # Apply the window function to the frames after post-processing
    # (take the real to ensure real values)
    audio_mdct = 2 * real(audio_mdct .* postprocessing_array) .* window_function

    # Loop over the time frames
    i = 0
    for j = 1:number_times

        # Recover the signal with the time-domain aliasing cancellation (TDAC) principle
        audio_signal[i+1:i+window_length] =
            audio_signal[i+1:i+window_length] + audio_mdct[:, j]
        i = i + step_length

    end

    # Remove the zero-padding at the start and at the end of the signal
    audio_signal = audio_signal[step_length+1:end-step_length]

end

"""
    window_function = hamming(window_length, window_sampling="symmetric")

Compute the Hamming window.

# Arguments:
- `window_length::Integer`: the window length in samples.
- `window_sampling::Char="symmetric"`: the window sampling method ("symmetric" or "periodic").
- `window_function::Float`: the window function (window_length,).
"""
function hamming(window_length, window_sampling = "symmetric")

    # Compute the Hamming window, symmetric for filter design and periodic for spectral analysis
    if window_sampling == "symmetric"
        window_function =
            0.54 .-
            0.46 * cos.(2 * pi * (0:window_length-1) / (window_length - 1))
    elseif window_sampling == "periodic"
        window_function =
            0.54 .- 0.46 * cos.(2 * pi * (0:window_length-1) / window_length)
    end

end

"""
    plot_object = zaf.sigplot(audio_signal, sampling_frequency, xtick_step=1)

Plot a signal in seconds.

# Arguments:
- `audio_signal::Float`: the audio signal (number_samples, number_channels).
- `sampling_frequency::Float`: the sampling frequency from the original signal in Hz.
- `xtick_step::Integer=1`: the step for the x-axis ticks in seconds (default: 1 second).
- `plot_object:Plots:` the plot object.
"""
function sigplot(audio_signal, sampling_frequency, xtick_step = 1)

    # Get the number of samples
    number_samples = size(audio_signal, 1)

    # Prepare the tick locations and labels for the x-axis
    xtick_locations = [
        xtick_step*sampling_frequency:xtick_step*sampling_frequency:number_samples;
    ]
    xtick_labels = convert(
        Array{Int},
        [xtick_step:xtick_step:number_samples/sampling_frequency;],
    )

    # Plot the signal in seconds
    plot_object = plot(
        audio_signal,
        legend = false,
        fmt = :png,
        xlims = (0, number_samples),
        xticks = (xtick_locations, xtick_labels),
        xlabel = "Time (s)",
    )

end

"""
    plot_object = zaf.specshow(audio_spectrogram, number_samples, sampling_frequency, xtick_step=1, ytick_step=1000)

Display a spectrogram in dB, seconds, and Hz.

# Arguments:
- `audio_spectrogram::Float`: the audio spectrogram (without DC and mirrored frequencies) (number_frequencies, number_times).
- `number_samples::Integer`: the number of samples from the original signal.
- `sampling_frequency::Float`: the sampling frequency from the original signal in Hz.
- `xtick_step::Integer=1`: the step for the x-axis ticks in seconds (default: 1 second).
- `ytick_step::Integer=1000`: the step for the y-axis ticks in Hz (default: 1000 Hz).
- `plot_object:Plots:` the plot object.
"""
function specshow(
    audio_spectrogram,
    number_samples,
    sampling_frequency,
    xtick_step = 1,
    ytick_step = 1000,
)

    # Get the number of frequency channels and time frames
    number_frequencies, number_times = size(audio_spectrogram)

    # Derive the number of Hertz and seconds
    number_hertz = sampling_frequency / 2
    number_seconds = number_samples / sampling_frequency

    # Derive the number of time frames per second and the number of frequency channels per Hz
    time_resolution = number_times / number_seconds
    frequency_resolution = number_frequencies / number_hertz

    # Prepare the tick locations and labels for the x-axis
    xtick_locations =
        [xtick_step*time_resolution:xtick_step*time_resolution:number_times;]
    xtick_labels = convert(Array{Int}, [xtick_step:xtick_step:number_seconds;])

    # Prepare the tick locations and labels for the y-axis
    ytick_locations = [
        ytick_step*frequency_resolution:ytick_step*frequency_resolution:number_frequencies;
    ]
    ytick_labels = convert(Array{Int}, [ytick_step:ytick_step:number_hertz;])

    # Display the spectrogram in dB, seconds, and Hz
    plot_object = heatmap(
        20 * log10.(audio_spectrogram),
        fillcolor = :jet,
        legend = false,
        fmt = :png,
        xticks = (xtick_locations, xtick_labels),
        yticks = (ytick_locations, ytick_labels),
        xlabel = "Time (s)",
        ylabel = "Frequency (Hz)",
    )

end

"""
    plot_object = zaf.melspecshow(mel_spectrogram, number_samples, sampling_frequency, window_length, xtick_step=1)

Display a mel spectrogram in dB, seconds, and Hz.

# Arguments:
- `mel_spectrogram::Float`: the mel spectrogram (number_mels, number_times).
- `number_samples::Integer`: the number of samples from the original signal.
- `sampling_frequency::Float`: the sampling frequency from the original signal in Hz.
- `window_length::Integer`: the window length from the Fourier analysis in number of samples.
- `xtick_step::Integer=1`: the step for the x-axis ticks in seconds (default: 1 second).
- `plot_object:Plots:` the plot object.
"""
function melspecshow(
    mel_spectrogram,
    number_samples,
    sampling_frequency,
    window_length,
    xtick_step = 1,
)

    # Get the number of mels and time frames
    number_mels, number_times = size(mel_spectrogram)

    # Derive the number of seconds and the number of time frames per second
    number_seconds = number_samples / sampling_frequency
    time_resolution = number_times / number_seconds
    
    # Derive the minimum and maximum mel
    minimum_mel = 2595 * log10(1 + (sampling_frequency / window_length) / 700)
    maximum_mel = 2595 * log10(1 + (sampling_frequency / 2) / 700)

    # Compute the mel scale (linearly spaced)
    mel_scale = range(minimum_mel, stop = maximum_mel, length = number_mels)

    # Derive the Hertz scale (log spaced)
    hertz_scale = 700 .* (10 .^ (mel_scale / 2595) .- 1)
    
    # Prepare the tick locations and labels for the x-axis
    xtick_locations =
        [xtick_step*time_resolution:xtick_step*time_resolution:number_times;]
    xtick_labels = convert(Array{Int}, [xtick_step:xtick_step:number_seconds;])

    # Prepare the tick locations and labels for the y-axis
    ytick_locations = [1:8:number_mels;]
    ytick_labels = convert(Array{Int}, round.(hertz_scale[1:8:number_mels]))

    # Display the mel spectrogram in dB, seconds, and Hz
    plot_object = heatmap(
        20 * log10.(mel_spectrogram),
        fillcolor = :jet,
        legend = false,
        fmt = :png,
        xticks = (xtick_locations, xtick_labels),
        yticks = (ytick_locations, ytick_labels),
        xlabel = "Time (s)",
        ylabel = "Frequency (Hz)",
    )

end

"""
    plot_object = zaf.mfccshow(audio_mfcc, number_samples, sampling_frequency, xtick_step=1)

Display MFCCs in seconds.

# Arguments:
- `audio_mfcc::Float`: the audio MFCC (number_coefficients, number_times).
- `number_samples::Integer`: the number of samples from the original signal.
- `sampling_frequency::Integer`: the sampling frequency from the original signal in Hz.
- `xtick_step::Integer=1`: the step for the x-axis ticks in seconds (default: 1 second).
- `plot_object:Plots:` the plot object.
"""
function mfccshow(audio_mfcc, number_samples, sampling_frequency, xtick_step = 1)

    # Get the number of time frames
    number_times = size(audio_mfcc, 2)

    # Derive the number of seconds and the number of time frames per second
    number_seconds = number_samples / sampling_frequency
    time_resolution = number_times / number_seconds

    # Prepare the tick locations and labels for the x-axis
    xtick_locations =
        [xtick_step*time_resolution:xtick_step*time_resolution:number_times;]
    xtick_labels = convert(
        Array{Int},
        [xtick_step:xtick_step:number_times/time_resolution;],
    )

    # Display the MFCCs in seconds
    plot_object = heatmap(
        audio_mfcc,
        fillcolor = :jet,
        legend = false,
        fmt = :png,
        xticks = (xtick_locations, xtick_labels),
        xlabel = "Time (s)",
        ylabel = "Coefficients",
    )

end

"""
    plot_object = zaf.cqtspecshow(cqt_spectrogram, time_resolution, frequency_resolution, minimum_frequency, xtick_step=1)

Display a CQT spectrogram in dB and seconds, and Hz.

# Arguments:
- `cqt_spectrogram::Float`: the CQT spectrogram (number_frequencies, number_times).
- `time_resolution::Integer`: number of time frames per second.
- `octave_resolution::Integer`: number of frequency channels per octave.
- `minimum_frequency::Float`: the minimum frequency in Hz.
- `xtick_step::Integer=1`: the step for the x-axis ticks in seconds (default: 1 second).
- `plot_object:Plots:` the plot object.
"""
function cqtspecshow(
    cqt_spectrogram,
    time_resolution,
    octave_resolution,
    minimum_frequency,
    xtick_step = 1,
)

    # Get the number of frequency channels and time frames
    number_frequencies, number_times = size(cqt_spectrogram)

    # Prepare the tick locations and labels for the x-axis
    xtick_locations =
        [xtick_step*time_resolution:xtick_step*time_resolution:number_times;]
    xtick_labels = convert(
        Array{Int},
        [xtick_step:xtick_step:number_times/time_resolution;],
    )

    # Prepare the tick locations and labels for the y-axis
    ytick_locations = [0:octave_resolution:number_frequencies;]
    ytick_labels = convert(
        Array{Int},
        minimum_frequency * 2 .^ (ytick_locations / octave_resolution),
    )

    # Display the CQT spectrogram in dB, seconds, and Hz
    plot_object = heatmap(
        20 * log10.(cqt_spectrogram),
        fillcolor = :jet,
        legend = false,
        fmt = :png,
        xticks = (xtick_locations, xtick_labels),
        yticks = (ytick_locations, ytick_labels),
        xlabel = "Time (s)",
        ylabel = "Frequency (Hz)",
    )

end

"""
    plot_object = zaf.cqtchromshow(cqt_chromagram, time_resolution, xtick_step=1)

Display a CQT chromagram in seconds.

# Arguments:
- `cqt_chromagram::Float`: the CQT chromagram (number_chromas, number_times).
- `time_resolution::Integer`: the number of time frames per second.
- `xtick_step::Integer=1`: the step for the x-axis ticks in seconds (default: 1 second).
- `plot_object:Plots:` the plot object.
"""
function cqtchromshow(cqt_chromagram, time_resolution, xtick_step = 1)

    # Get the number of time frames
    number_times = size(cqt_chromagram, 2)

    # Prepare the tick locations and labels for the x-axis
    xtick_locations =
        [xtick_step*time_resolution:xtick_step*time_resolution:number_times;]
    xtick_labels = convert(
        Array{Int},
        [xtick_step:xtick_step:number_times/time_resolution;],
    )

    # Display the CQT chromagram in seconds
    plot_object = heatmap(
        cqt_chromagram,
        fillcolor = :jet,
        legend = false,
        fmt = :png,
        xticks = (xtick_locations, xtick_labels),
        xlabel = "Time (s)",
        ylabel = "Chroma",
    )

end

end
