clc
clear all
close all
% Assume 0 degrees due north, with positive clockwise orientation
% Positive phase difference indicates that antennna 1 recieves the waveform
% before the respective antenna (2,3, or 4)
% 
% Fix antenna 1 at the origin. From the x-axis, rotate antenna 2 60 degrees,
% rotate antenna 3 90 degrees, and rotate antenna 4 120 degrees. The
% optimal antenna spacing is still to be determined.
% 
%                         | 0 deg         --> +
%                         |
%                         |
%                         |
%             ---------- A1 ----------
%                   -60 ' | `  60
%                      '  |  `
%                     '   A3  `
%                  A4'         `A2

%%

num_antennas = 10;
num_sig_sources = 2;
RF_freq = 600e6;        % freq of incoming signal
lambda = 3e8/RF_freq;   % wavelength of incoming signal
fs = 100e6;             % sampling frequency
samples = 2000;          % number of snapshots
cos_freq = 1e6;         % freq to downconvert to

ant_spacings = [.25, .98, .91]; % antenna spacing

% generate_lookup_table just calculates a set of theoretical phase
% differences corresponding to some random AoA for the purpose of signal
% generation.
n = randi([0,360],1);
[AoAs, PDs1to2, PDs1to3, PDs1to4] = generate_lookup_table(-90, 90, .5, lambda, ant_spacings(1), ant_spacings(2), ant_spacings(3));
AoA = AoAs(n);

% generate signal
t = linspace(0,samples/fs,samples);
% recived_signal_ant_1 = cos(2*pi.*t.*cos_freq);
% recived_signal_ant_2 = cos(2*pi.*t.*cos_freq + PDs1to2(n));
% recived_signal_ant_3 = cos(2*pi.*t.*cos_freq + PDs1to3(n));
% recived_signal_ant_4 = cos(2*pi.*t.*cos_freq + PDs1to4(n));

% Change signals to be complex
recived_signal_ant_1 = exp(1i*2*pi.*t.*cos_freq);
recived_signal_ant_2 = exp(1i*(2*pi.*t.*cos_freq + PDs1to2(n)));
recived_signal_ant_3 = exp(1i*(2*pi.*t.*cos_freq + PDs1to3(n)));
recived_signal_ant_4 = exp(1i*(2*pi.*t.*cos_freq + PDs1to4(n)));


% build sample correlation matrix
cov_mat = zeros(4,4);
for i = 1:samples
    X = [recived_signal_ant_1(i);
         recived_signal_ant_2(i); 
         recived_signal_ant_3(i);
         recived_signal_ant_4(i)];
    cov_mat = cov_mat + X*X';
end

cov_mat = cov_mat*(1/samples);

% grab eigenvalues and eigenvectors
% [eig_vecs, eig_vals] = eig(cov_mat);

% generate matrix w/ noise eigenvector columns. An eigenvector corresponds
% to noise if the associated eigenvalue is close to 0. 
% col_num = 1;
% for i = 1:4
%     if eig_vals(i,i) < 1e-4
%         noise_vecs(:,col_num) = eig_vecs(:,i);
%         col_num = col_num + 1;
%     end
% end

% Use svd() instead
[eig_vecs, eig_vals] = svd(cov_mat);
noise_vecs = eig_vecs(:,num_sig_sources+1:end);

% walk t
theta = linspace(-90,90,200);
spectrum_func = zeros(1,100);
for i = 1:length(theta)

    % steering vector - only taking real part as of now (Re(e^i*phasediff) = cos(phase_diff))
    steering_vec = [1;
%         cos((-2*pi*ant_spacings(1)/lambda)*sin(theta(i)*(pi/180) - pi/3));
%         cos((-2*pi*ant_spacings(2)/lambda)*sin(theta(i)*(pi/180) - pi/2));
%         cos((-2*pi*ant_spacings(3)/lambda)*sin(theta(i)*(pi/180) - 2*pi/3))];

    % Make the steering vectors complex
    exp((1i*(-2*pi*ant_spacings(1)/lambda)*sin(theta(i)*(pi/180) - pi/3)));
        exp((1i*(-2*pi*ant_spacings(2)/lambda)*sin(theta(i)*(pi/180) - pi/2)));
        exp((1i*(-2*pi*ant_spacings(3)/lambda)*sin(theta(i)*(pi/180) - 2*pi/3)))];
    spectrum_func(i) = 1/((steering_vec'*noise_vecs)*noise_vecs'*steering_vec);
end

% grab theta corresponding to highest peak, round to nearest half integer
[mx,argmax] = max(spectrum_func);
estimated_AoA = floor(2*theta(argmax)+.5)/2;
fprintf("Estimated AoA: %.1f \t Hard-coded AoA: %.1f.\n", estimated_AoA, AoA)

plot(theta,20*log10(abs((spectrum_func))))
title('MUSIC Psuedospectrum')
xlabel('Angle of Arrival')





function [AoAs, PDs1to2, PDs1to3, PDs1to4] = generate_lookup_table(start, stop, step, lambda, d1to2, d1to3, d1to4)
% d1to2 = spacing between antennas 1 and 2

num_pts = (stop-start)/step + 1;
AoAs = zeros(1,num_pts);
PDs1to2 = zeros(1,num_pts);
PDs1to3 = zeros(1,num_pts);
PDs1to4 = zeros(1,num_pts);

for i = 1:num_pts 
    AoAs(i) = start + (i-1)*step;

    % phase differences
    value_to_push_to_2 = (-2*pi*d1to2/lambda)*sin(AoAs(i)*(pi/180) - pi/3);
    value_to_push_to_3 = (-2*pi*d1to3/lambda)*sin(AoAs(i)*(pi/180) - pi/2);
    value_to_push_to_4 = (-2*pi*d1to4/lambda)*sin(AoAs(i)*(pi/180) - 2*pi/3);

    % increment/decrement until in -180 to 180
    if value_to_push_to_2 > 2*pi
        n = ceil(value_to_push_to_2/(2*pi) - 1);
        value_to_push_to_2 = value_to_push_to_2 - 2*pi*n;
    end

    if value_to_push_to_2 > pi
        value_to_push_to_2 = value_to_push_to_2 - 2*pi;
    end

    if value_to_push_to_3 > 2*pi
        n = ceil(value_to_push_to_3/(2*pi) - 1);
        value_to_push_to_3 = value_to_push_to_3 - 2*pi*n;
    end

    if value_to_push_to_3 > pi
        value_to_push_to_3 = value_to_push_to_3 - 2*pi;
    end

    if value_to_push_to_4 > 2*pi
        n = ceil(value_to_push_to_4/(2*pi) - 1);
        value_to_push_to_4 = value_to_push_to_4 - 2*pi*n;
    end

    if value_to_push_to_4 > pi
        value_to_push_to_4 = value_to_push_to_4 - 2*pi;
    end

    if value_to_push_to_2 < -2*pi
        n = ceil(-1 - value_to_push_to_2/(2*pi));
        value_to_push_to_2 = value_to_push_to_2 + 2*pi*n;
    end

    if value_to_push_to_2 < -pi
        value_to_push_to_2 = value_to_push_to_2 + 2*pi;
    end

    if value_to_push_to_3 < -2*pi
        n = ceil(-1 - value_to_push_to_3/(2*pi));
        value_to_push_to_3 = value_to_push_to_3 + 2*pi*n;
    end

    if value_to_push_to_3 < -pi
        value_to_push_to_3 = value_to_push_to_3 + 2*pi;
    end
    if value_to_push_to_4 < -2*pi
        n = ceil(-1 - value_to_push_to_4/(2*pi));
        value_to_push_to_4 = value_to_push_to_4 + 2*pi*n;
    end

    if value_to_push_to_4 < -pi
        value_to_push_to_4 = value_to_push_to_4 + 2*pi;
    end

    if value_to_push_to_4 < -2*pi
        n = ceil(-1 - value_to_push_to_4/(2*pi));
        value_to_push_to_4 = value_to_push_to_4 + 2*pi*n;
    end

    if value_to_push_to_4 < -pi
        value_to_push_to_4 = value_to_push_to_4 + 2*pi;
    end

    PDs1to2(i) = value_to_push_to_2;
    PDs1to3(i) = value_to_push_to_3;
    PDs1to4(i) = value_to_push_to_4;
end
end