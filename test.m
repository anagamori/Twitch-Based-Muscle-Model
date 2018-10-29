close all
clear all
clc

Fs = 1000; %sampling frequency
time = 0:1/Fs:5; %simulation time

T_1 = 0.01;
T_2 = 0.03;

spike = zeros(1,length(time));
x_1_vec = zeros(1,length(time));
x_2_vec = zeros(1,length(time));
spike(1*Fs) = 1;

x_1 = 0;
x_2 = 0;
for t = 1:length(time)
   x_1_dot = (spike(t)-x_1)/T_1;
   x_1 = x_1_dot/Fs + x_1;
   x_2_dot = (spike(t)-x_2)/T_2;
   x_2 = x_2_dot/Fs + x_2;
   
   x_1_vec(t) = x_1;
   x_2_vec(t) = x_2;
   
   f = x_2/(x_2+x_1);
   f_vec(t) = f;
end

[pks,locs_peak] = max(x_2_vec);
CT = locs_peak-1*Fs;

peak_half = pks/2;
[~,HRT] = min(abs(x_2_vec(locs_peak:end)-peak_half));
locs_hrt = locs_peak+HRT;

figure(1)
subplot(2,1,1)
plot(time,x_1_vec)
subplot(2,1,2)
plot(time,x_2_vec)
hold on 
plot([time(locs_peak) time(locs_peak)],[0 1],'--')
text(time(locs_peak),pks,num2str(CT))
hold on
plot([time(locs_hrt) time(locs_hrt)],[0 1],'--')
text(time(locs_hrt),pks,num2str(HRT))

figure(2)
plot(time,f_vec)