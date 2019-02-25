%==========================================================================
% new_model_test_f2a.m
% Author: Akira Nagamori
% Last update: 10/31/18
%
%==========================================================================

%close all
clear all
clc

%==========================================================================
% Simulation parameters
FR_test = [2 5 10 15 20 30 50 100]; %10:10:100];
Fs = 1000; %sampling frequency
T = 1/Fs;
time = 0:1/Fs:4; %simulation time

simulation_condition = 2;

%==========================================================================
% initialization
mean_exc = zeros(1,length(FR_test));
p2p_exc = zeros(1,length(FR_test));

for f = 1:length(FR_test)
    FR = FR_test(f);
    %==========================================================================
    % Generate spike train
    spike = zeros(1,length(time));
    
    % Generate spike train
    temp =spikeTrainGenerator(0:1/Fs:2,Fs,FR);
    spike(1*Fs:3*Fs) = temp;
    
    %==========================================================================
    % Model parameters
    T_1 = 50*Fs/1000;
    T_2 = 30*Fs/1000;
    T_3 = 50*Fs/1000;
    T_4_0 = 100*Fs/1000;
    
    tau_1 = 0.003*Fs/1000;
    tau_2 = 0.024*Fs/1000;    
    
    n = 3;
    k = 0.3;
    N = 2;
    
    %======================================================================
    % initialization
    R = zeros(1,length(time));
    x_1_int_vec = zeros(1,length(time));
    x_1_vec = zeros(1,length(time));
    y_int_vec = zeros(1,length(time));
    y_vec = zeros(1,length(time));
    x_1 = 0;
    y_int = 0;
    y = 0;
    
    
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
        
        x_1_dot = R(t)^N*T_1-x_1^N*T_2;
        x_1 = x_1_dot/Fs + x_1;
        %T_4 = T_4_0*(1-0.01*y)^2;
        T_4 = T_4_0/(1+5*y);
        %T_4 = T_4_0;
        y_dot = x_1^N*T_3-y_int^N*T_4;
        y_int = y_dot/Fs + y_int;
        y = y_int^n/(y_int^n+k^n);
        %R_vec(t) = R;
        x_1_vec(t) = x_1;
        y_int_vec(t) = y_int;
        y_vec(t) = y;
        
    end
    
    mean_exc(f) = mean(y_vec(2*Fs:3*Fs));
    p2p_exc(f) = max(y_vec(2*Fs:3*Fs))-min(y_vec(2*Fs:3*Fs));
    
    
    figure(1)
    plot(time,y_vec)
    hold on
    
end

fusion = 1-p2p_exc/p2p_exc(1);
figure(2)
plot(FR_test,mean_exc,'LineWidth',1)
hold on

figure(3)
plot(FR_test,fusion,'LineWidth',1)
hold on

figure(4)
plot(mean_exc,fusion,'LineWidth',1)
hold on

function spikeTrain = spikeTrainGenerator(t,Fs,freq)

spikeTrain = zeros(1,length(t));
ISI = round(1/freq*Fs);
numSpikes = round(length(t)/ISI);
index = [1:numSpikes]*ISI;
index(index>length(t)) = [];
spikeTrain(index) = 1;
spikeTrain(1) = 1;

end
