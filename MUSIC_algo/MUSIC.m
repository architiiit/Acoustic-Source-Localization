%% MUSIC Algorithm
%% General
clc
close all
clear all

N=200;%Samples
doa=[30 -60]/180*pi;%Angles At which sources are placed
w=[pi/4 pi/3]';%frequency
M=10;%Array Numbers
lambda=150;%spacing between array elements
d=lambda/2;%array element space
snr=20;%Signal to noise ratio

P=length(w);
B=zeros(P,M);

for k=1:1:P
 B(k,:)=exp(-j*2*pi*d*sin(doa(k))/lambda*[0:M-1]);
end


B=B';
xx=2*exp(j*(w*[1:N]));
x=B*xx;
x=x+awgn(x,snr);%Gaussin noise
R=x*x';
[U,V]=eig(R);
UU=U(:,1:M-P);%noise sub space
theta=-90:0.5:90;
for ii=1:length(theta)
 AA=zeros(1,length(M));
 for jj=0:M-1
 AA(1+jj)=exp(-j*2*jj*pi*d*sin(theta(ii)/180*pi)/lambda);
 end
 WW=AA*UU*UU'*AA';
 Pmusic(ii)=abs(1/WW);
end
Pmusic=10*log10(Pmusic/max(Pmusic));%spatial spectrum
plot(theta,Pmusic,'-k','linewidth',2.0)
xlabel('DOA \theta/degree')
ylabel('Power Spectrum/dB')
title('MUSIC algorithm for multiple sources')
grid on