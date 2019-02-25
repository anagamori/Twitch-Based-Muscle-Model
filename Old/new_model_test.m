%==========================================================================
% new_model_test.m
% Author: Akira Nagamori
% Last update: 11/1/18
%==========================================================================

close all
clear all
clc

%==========================================================================
% Simulation parameters
FR = 1;
Fs = 1000; %sampling frequency
T = 1/Fs;
time = 0:1/Fs:4; %simulation time

simulation_condition = 2;
%==========================================================================
% Generate spike train
spike = zeros(1,length(time));
if simulation_condition == 1
    % Generate a pulse
    spike(1*Fs) = 1;
else
    % Generate spike train
    temp =spikeTrainGenerator(0:1/Fs:2,Fs,FR);
    spike(1*Fs:3*Fs) = temp;
end

%==========================================================================
% Model parameters
T_1 = 50*Fs/1000;
T_2_0 = 50*Fs/1000;
T_3 = 1000*Fs/1000;
T_4_0 = 100*Fs/1000;

tau_1 = 0.003*Fs/1000;
tau_2 = 0.024*Fs/1000;

n = 3;
k = 0.3;
N = 2;
%==========================================================================
% initialization
R = zeros(1,length(time));
x_1_int_vec = zeros(1,length(time));
x_1_vec = zeros(1,length(time));
y_int_vec = zeros(1,length(time));
y_vec = zeros(1,length(time));
x_1 = 0.05;
y_int = 0;
y = 0.0026;

for t = 1:length(time)
    spike_temp = zeros(1,length(time));
    if spike(t) == 1
        spike_temp(t) = 1;
        R_temp_1 = conv(spike_temp,1-exp(-time/tau_1));
        R_temp_1 = R_temp_1(1:length(time));
        R_temp_2 = conv(spike_temp,exp(-time/tau_2));
        R_temp_2 = R_temp_2(1:length(time));
        R = R_temp_1.*R_temp_2+R;
        R = R(1:length(time));
    end
    
    T_2 = T_2_0/(1+5*y);
    x_1_dot = R(t)^N*T_1-x_1^N*T_2;
    x_1 = x_1_dot/Fs + x_1;
    %T_4 = T_4_0*(1-0.01*y)^2;
    T_4 = T_4_0*x_1^2;
    %T_4 = T_4_0;
    y_dot = x_1^N*T_3-y_int^N*T_4;
    y_int = y_dot/Fs + y_int;
    y = y_int^n/(y_int^n+k^n);
    %R_vec(t) = R; 
    x_1_vec(t) = x_1;
    y_int_vec(t) = y_int;
    y_vec(t) = y;
    
end

if simulation_condition == 1
    % Twitch analysis
    [pks,locs_peak] = max(y_vec);
    pks
    CT = (locs_peak-1*Fs)*1000/Fs
    
    peak_half = pks/2;
    [~,HRT] = min(abs(y_vec(locs_peak:end)-peak_half));
    locs_hrt = locs_peak+HRT;
    HRT = HRT*1000/Fs
else
    mean_exc = mean(y_vec(1*Fs:2*Fs))
end


figure(1)
plot(time,y_vec)
hold on
% plot(time,y_int_vec)

if simulation_condition == 1
    hold on
    plot([time(locs_peak) time(locs_peak)],[0 1],'--')
    text(time(locs_peak),pks,num2str(CT))
    hold on
    plot([time(locs_hrt) time(locs_hrt)],[0 1],'--')
    text(time(locs_hrt),peak_half,num2str(HRT))
    text(time(locs_peak)-time(locs_peak)*0.3,pks,num2str(pks))
end

function spikeTrain = spikeTrainGenerator(t,Fs,freq)

spikeTrain = zeros(1,length(t));
ISI = round(1/freq*Fs);
numSpikes = round(length(t)/ISI);
index = [1:numSpikes]*ISI;
index(index>length(t)) = [];
spikeTrain(index) = 1;
spikeTrain(1) = 1;

end
