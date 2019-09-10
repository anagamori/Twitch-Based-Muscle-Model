close all
clear all
clc

%==========================================================================
% Simulation parameters
Fs = 1000; %sampling frequency
T = 1/Fs;
time = 0:1/Fs:5; %simulation time

tau_1 = 0.005;


FR = 100;
spike = zeros(1,length(time));
temp =spikeTrainGenerator(0:1/Fs:3,Fs,FR);
spike(1*Fs:4*Fs) = temp;

spike_temp = zeros(1,length(time));
R_temp = exp(-time/tau_1);
R = zeros(1,length(time));
R_2_vec = zeros(1,length(time));
R_2 = 0;

alpha = 3;
A = 3;

for t = 1:length(time)
    
    spike_temp = zeros(1,length(time));
    if spike(t) == 1
        spike_temp(t) = 1;
        temp = conv(spike_temp,R_temp);
        R = R + temp(1:length(time));
    end
    
    R_2 = spike(t) + exp(-T/tau_1)*R_2;
    R_2_vec(t) = R_2;
end

figure(1)
plot(time,R)
hold on
plot(time,R_2_vec,'--')

function spikeTrain = spikeTrainGenerator(t,Fs,freq)

spikeTrain = zeros(1,length(t));
ISI = round(1/freq*Fs);
numSpikes = round(length(t)/ISI);
index = [1:numSpikes]*ISI;
index(index>length(t)) = [];
spikeTrain(index) = 1;
spikeTrain(1) = 1;

end