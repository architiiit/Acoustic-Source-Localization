% Generate synthetic data
% Example inputs
clc
clear all
close all

% YY = randn(8, 1000);  % Example received signal snapshots (8 sensors, 1000 samples)
% k = 3;                % Desired number of DOAs
% B = 1;                % Number of blocks
% 
% % Perform Root MUSIC algorithm
% [sr, alr] = rootMUSIC(YY, k, B);
% 
% % Display results
% disp('Estimated DOAs (in radians):');
% disp(sr);
% disp('All roots of the polynomial (in radians):');
% disp(alr);
% % Plot estimated DOAs
% figure;
% subplot(2,1,1);
% stem(sr, ones(size(sr)), 'o', 'LineWidth', 2);
% title('Estimated DOAs');
% xlabel('Angle (radians)');
% ylabel('Magnitude');
% grid on;
% 
% % Plot all roots of the polynomial
% subplot(2,1,2);
% stem(alr, ones(size(alr)), 'x', 'LineWidth', 2);
% title('All Roots of the Polynomial');
% xlabel('Angle (radians)');
% ylabel('Magnitude');
% grid on;
% Parameters
% Parameters
N = 128; % Number of samples
M = 8; % Number of sensors
d = 0.5; % Half wavelength as array element spacing
f_center = 1000; % Center frequency in Hz
SNRs = [10, 20, 30]; % Signal-to-noise ratios in dB
incident_angles = [60, 30, 75]; % Incident angles in degrees

% Generate synthetic signals from narrowband sources
signal_sources = zeros(M, N);
for i = 1:length(incident_angles)
    theta_rad = deg2rad(incident_angles(i));
    s_i = exp(1j*2*pi*f_center*d*sin(theta_rad)*(0:M-1)') / sqrt(M);
    signal_sources = signal_sources + repmat(s_i, 1, N);
end

% Perform simulation for different SNRs
estimated_DOAs = cell(length(SNRs), 1);
all_roots = cell(length(SNRs), 1);

for snr_idx = 1:length(SNRs)
    SNR = SNRs(snr_idx);
    noise_power = 10^(-SNR/10);
    noise = sqrt(noise_power/2) * (randn(M, N) + 1j*randn(M, N));
    received_data = signal_sources + noise;
    
    % Perform Root MUSIC algorithm
    [sr, alr] = rootMUSIC(received_data, length(incident_angles), 1);
    
    estimated_DOAs{snr_idx} = sr;
    all_roots{snr_idx} = alr;
end

% Plot results
figure;
for snr_idx = 1:length(SNRs)
    subplot(length(SNRs), 2, 2*snr_idx - 1);
    stem(rad2deg(estimated_DOAs{snr_idx}), ones(size(estimated_DOAs{snr_idx})), 'o', 'LineWidth', 2);
    title(['Estimated DOAs (SNR = ' num2str(SNRs(snr_idx)) ' dB)']);
    xlabel('Angle (degrees)');
    ylabel('Magnitude');
    grid on;
    
    subplot(length(SNRs), 2, 2*snr_idx);
    stem(rad2deg(all_roots{snr_idx}), ones(size(all_roots{snr_idx})), 'x', 'LineWidth', 2);
    title(['All Roots of the Polynomial (SNR = ' num2str(SNRs(snr_idx)) ' dB)']);
    xlabel('Angle (degrees)');
    ylabel('Magnitude');
    grid on;
end



function [sr, alr] = rootMUSIC(YY, k, B)
% root_MUSIC: Root MUSIC algorithm for direction-of-arrival (DOA) estimation
% 
% Inputs:
%   YY: Signal Snapshots/Samples
%   k: Number of true features (desired number of DOAs)
%   B: Number of blocks
%
% Outputs:
%   sr: Estimated DOAs (in radians)
%   alr: All roots of the polynomial (in radians)
%
% Note: This function assumes that the input YY is a matrix where each column
% represents a snapshot or sample of received signal data.

[L, ~] = size(YY);
R = YY * YY' / (L * B); % Covariance of data matrix
[eigvec, ~] = eig(R);

%% Computing coefficients for polynomial construction
En = eigvec(:, 1:L - k); % noise eigenspace
A = En * En';
cf = zeros(L - 1, 1);
for i = 1:L - 1
    cf(i) = sum(diag(A, i));
end
% Complete set of coefficients
CF = [flipud(cf); sum(diag(A)); conj(cf)];

% Computing roots
rts = roots(CF);
% Finding roots inside the unit circle and refining the roots
candidates = abs(rts) <= 1;
true_rts = rts(candidates);
[~, ind] = maxk(abs(true_rts), k);
sr = angle(true_rts(ind));

% Convert estimated DOAs to the range [0, 2*pi)
sr = wrapTo2Pi(sr);

% Compute angles of all roots
alr = wrapTo2Pi(angle(rts));
end

