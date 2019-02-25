%--------------------------------------------------------------------------
% twitch_analysis.m
% Author: Akira Nagamori
% Last update: 10/16/18
% Descriptions: analyze effects of force-dependent time constant on contraction time and
% half relaxation time 
%--------------------------------------------------------------------------
%close all
clc

Fs = 1000; %sampling frequency
time = 0:1/Fs:5; %simulation time
stim = zeros(1,length(time)); %stimulus train
stim(1*Fs) = 1; %single pulse stimulus

a = 10;
b = 5;

h = exp(-a*time); %first linear filter
g = exp(-b*time); %second linear filter whose time constant depends on the force level (Bobet & Stein, 1998)

q = conv(stim,h); %convolution
f = conv(q,g); %convolution
x = conv(f,g);

[pks,locs_peak] = max(x);
CT = locs_peak-1*Fs;

peak_half = pks/2;
[~,HRT] = min(abs(x(locs_peak:end)-peak_half));
locs_hrt = locs_peak+HRT;

figure()
plot(time,x(1:length(time)),'LineWidth',1)
hold on 
plot([time(locs_peak) time(locs_peak)],[0 pks+10],'--')
text(time(locs_peak),pks,num2str(CT))
hold on
plot([time(locs_hrt) time(locs_hrt)],[0 pks+10],'--')
text(time(locs_hrt),pks,num2str(HRT))






