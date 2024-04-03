% FFT based DOA estimation
%% single source, no window function
% clc
% close all
% clear all
% 
% lamda=1;
% d=lamda/2;
% N=8; % array number
% M=1; % number of sources
% rad=pi/180;
% 
% theta_1=30.1; % DOA of source 1
% theta_1=theta_1*rad;
% phi_1=2*pi*(d*cos(theta_1)/lamda);
% 
% n=0:1:(N-1);
% n=n';
% 
% s_1=exp(1i*n*phi_1);
% SNR=20;
% 
% Ps=s_1'*s_1;
% Pn=Ps*10^(-SNR/10);
% 
% sigma=sqrt(Pn/(2*N));
% x=s_1+sigma*(randn(N,1)+1j*rand(N,1));
% 
% figure(1)
% subplot(2,1,1)
% plot(abs(s_1))
% 
% subplot(2,1,2)
% plot(abs(x))
% 
% z=1024;
% L=z*N;
% X=fft(x,L);
% Y=fftshift(abs(X).^2);
% phi_axis=-180:360/L:180-360/L;
% theta_axis=acosd(phi_axis/360*lamda/d);
% [A,m]=max(Y);
% 
% 
% figure(2)
% plot(phi_axis,10*log10(fftshift(Y)),'LineWidth',2);
% axis([-180 180 -20 1.1*max(10*log10(fftshift(Y)))])
% xlabel('DOA \theta/radians')
% ylabel('Power Spectrum/dB')
% title('Single Source DOA Estimation')
% grid on
% phi_0=phi_axis(m)
% 
% 
% figure(3)
% plot(theta_axis,10*log10(fftshift(Y)),'LineWidth',2);
% axis([0 180 -20 1.1*max(10*log10(fftshift(Y)))])
% xlabel('DOA \theta/degree')
% ylabel('Power Spectrum/dB')
% title('Single Source DOA Estimation')
% grid on
% theta_0=theta_axis(m)


% X_axis=acos(x_axis*rad/(2*pi)*lamda/d)/rad
% figure(2)
% plot(spec)
% axis([0 360 min(spec) 1.2*max(spec)])
% grid on


% Multiple targets
clc
close all
clear all
lamda=1;
d=lamda/2;
M=10; % array number
Source=[1 ; exp(1i*pi/4)]; % number of sources
rad=pi/180;
theta_1=30; % DOA of source 1
theta_1=theta_1*rad;
phi_1=2*pi*(d*cos(theta_1)/lamda);

theta_2=46.5; % DOA of source 2
theta_2=theta_2*rad;
phi_2=2*pi*(d*cos(theta_2)/lamda);

M_array=0:1:(M-1);
M_array=M_array';

s_1=exp(1i*M_array*phi_1);
s_2=exp(1i*M_array*phi_2)
;
s=[s_1 s_2];
ss=s*Source;
SNR=20;

x=ss+awgn(ss,SNR);

z=1024;
L=z*M;
X=fft(x,L);
Y=fftshift(abs(X).^2);
phi_axis=-180:360/L:180-360/L;
theta_axis=acosd(phi_axis/360*lamda/d);
[A,m]=max(Y);

figure(2)
xlabel('DOA \theta/degree')
ylabel('Power Spectrum/dB')
title('Multiple Sources DOA Estimation')
grid on
plot(phi_axis,10*log10(fftshift(Y)));
axis([-180 180 min(10*log10(fftshift(Y)))
1.1*max(10*log10(fftshift(Y)))])

phi_0=phi_axis(m)

figure(3)
plot(theta_axis,10*log10(fftshift(Y)));
axis([0 180 min(10*log10(fftshift(Y)))
1.1*max(10*log10(fftshift(Y)))])
xlabel('DOA \theta/degree')
ylabel('Power Spectrum/dB')
title('Multiple Sources DOA Estimation')
grid on
theta_0=theta_axis(m)