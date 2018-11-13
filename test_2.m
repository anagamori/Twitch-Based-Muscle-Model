close all
clear all
clc

Fs = 1000; %sampling frequency
time = 0:1/Fs:3; %simulation time

T_1 = 0.003;
T_2 = 0.025;

spike = zeros(1,length(time));
f = 1-exp(-time./T_1);
g = exp(-time/T_2);
x_1_vec = zeros(1,length(time));
x_2_vec = zeros(1,length(time));
spike(1*Fs) = 1;

x_1 = 0;
x_2 = 0;

x = conv(spike,f);
y = conv(spike,g);
z = x.*y;

z = z(1:length(time));

[pks,locs_peak] = max(z);
CT = locs_peak-1*Fs;

peak_half = pks/2;
[~,HRT] = min(abs(z(locs_peak:end)-peak_half));
locs_hrt = locs_peak+HRT;
% 
figure(1)
plot(time,z)
hold on 
plot([time(locs_peak) time(locs_peak)],[0 1],'--')
text(time(locs_peak),pks,num2str(CT))
hold on
plot([time(locs_hrt) time(locs_hrt)],[0 1],'--')
text(time(locs_hrt),pks,num2str(HRT))
% 
% figure(2)
% plot(time,f_vec)