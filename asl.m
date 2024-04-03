%code to perform mvdr and das
clc
clear all
close all

% Load the audio signal
[audio_signal, fs] = audioread('sample.wav');
M = 64; % Number of microphones
c = 1500; % Speed of sound in m/s
d = 0.5; % Distance between microphones in meters
t_delays = (0:M-1) * d / c;
angles_deg = rad2deg(asin(t_delays * c));
% Preallocate matrix for delayed signals
% delayed_signals = zeros(length(audio_signal), M);
% % Calculate delayed signals for DAS beamforming
% for mic = 1:M
% delayed_signals(:, mic) = circshift(audio_signal, round(t_delays(mic) * fs));
% end
% % Beamforming for DAS
% beamformed_signal_das = sum(delayed_signals, 2);
% beamformed_signal_das = beamformed_signal_das / max(abs(beamformed_signal_das));
% % Calculate noise power for DAS
% noise_power_das = var(audio_signal - beamformed_signal_das);
% beamformed_power_das = var(beamformed_signal_das);
% SNR_delay_sum = 10 * log10(beamformed_power_das / noise_power_das);
% % MVDR Beamformer
% Rxx = audio_signal' * audio_signal / length(audio_signal);
% epsilon = 1e-6;
% Rxx_inv = inv(Rxx + epsilon * eye(M));
% mvdr_weights = (Rxx_inv * ones(M, 1)) / (ones(1, M) * Rxx_inv * ones(M, 1));
% mvdr_signal = audio_signal * mvdr_weights';
% mvdr_signal = mvdr_signal / max(abs(mvdr_signal));
% % Calculate noise power for MVDR
% noise_power_mvdr = var(audio_signal - mvdr_signal);
% mvdr_power = var(mvdr_signal);
% SNR_mvdr = 10 * log10(mvdr_power / noise_power_mvdr);
% % Plotting and displaying SNR values
% figure;
% subplot(2,1,1);
% t_audio = (0:length(audio_signal)-1) / fs; % Time vector for original audio signal
% plot(t_audio, beamformed_signal_das);
% title('Beamformed Signal (DAS)');
% xlabel('Time (s)');
% ylabel('Amplitude');
% subplot(2,1,2);
% plot(t_audio, mvdr_signal);
% title('MVDR Signal');
% xlabel('Time (s)');
% ylabel('Amplitude');
% fprintf('SNR (Delay-and-Sum): %.2f dB\n', SNR_delay_sum);
% fprintf('SNR (MVDR): %.2f dB\n', SNR_mvdr);

% 
% 
% code to plot amplitude spectrum of beamformed signal
% 
% 
% % Load the audio signal
% [audio_signal, fs] = audioread('sample.wav');
% M = 64; % Number of microphones
% c = 1500; % Speed of sound in m/s
% d = 0.5; % Distance between microphones in meters
% t_delays = (0:M-1) * d / c;
% angles_deg = rad2deg(asin(t_delays * c));
% % Preallocate matrix for delayed signals
% delayed_signals = zeros(length(audio_signal), M);
% % Calculate delayed signals for DAS beamforming
% for mic = 1:M
% delayed_signals(:, mic) = circshift(audio_signal, round(t_delays(mic) * fs));
% end
% % Beamforming for DAS
% beamformed_signal_das = sum(delayed_signals, 2);
% beamformed_signal_das = beamformed_signal_das / max(abs(beamformed_signal_das));
% % MVDR Beamformer
% Rxx = audio_signal' * audio_signal / length(audio_signal);
% epsilon = 1e-6;
% Rxx_inv = inv(Rxx + epsilon * eye(M));
% mvdr_weights = (Rxx_inv * ones(M, 1)) / (ones(1, M) * Rxx_inv * ones(M, 1));
% mvdr_signal = audio_signal * mvdr_weights';
% mvdr_signal = mvdr_signal / max(abs(mvdr_signal));
% % Plotting the amplitude spectrum of beamformed signals
% figure;
% subplot(2,1,1);
% NFFT = 2^nextpow2(length(beamformed_signal_das));
% frequencies = fs / 2 * linspace(0, 1, NFFT/2 + 1);
% beamformed_das_fft = fft(beamformed_signal_das, NFFT);
% amplitude_spectrum_das = 2 * abs(beamformed_das_fft(1:NFFT/2 + 1));
% plot(frequencies, 20*log10(amplitude_spectrum_das));
% title('Amplitude Spectrum of Beamformed Signal (DAS)');
% xlabel('Frequency (Hz)');
% ylabel('Amplitude (dB)');
% grid on;
% subplot(2,1,2);
% mvdr_signal_fft = fft(mvdr_signal, NFFT);
% amplitude_spectrum_mvdr = 2 * abs(mvdr_signal_fft(1:NFFT/2 + 1));
% plot(frequencies, 20*log10(amplitude_spectrum_mvdr));
% title('Amplitude Spectrum of MVDR Signal');
% xlabel('Frequency (Hz)');
% ylabel('Amplitude (dB)');
% grid on;