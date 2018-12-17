%==========================================================================
% new_model_v2_test.m
% Author: Akira Nagamori
% Last update: 11/13/18
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
    x = 0; % proportion of Ca2+ bound to troponin (0-1)
    x_int = 0;
    y = 0; % porportion of cross-bridge sites available for binding (0-1)
    z = 0; % porportion of bidning sites that formed cross-bridges (0-1)
    
    T_x_1 = 400*Fs/1000; % forward rate constant for Ca2+ binding to troponin
    T_x_2 = 50*Fs/1000; % backward rate constant for Ca2+ release from troponin
    T_y_1 = 200*Fs/1000; % forward rate constant for binding sites becoming available
    T_y_2 = 80*Fs/1000; % backward rate constant for binding sites becoming unavailable
    T_z_1 = 10*Fs/1000; % forward rate constant for binding sites forming cross-bridge
    T_z_2 = 40*Fs/1000; % backward rate constant for binding sites detachecing
    
    K = 300;
    alpha = 10;
    
    n = 3;
    k = 0.3;
    N = 3;
    %==========================================================================
    % initialization
    
    
    x_vec = zeros(1,length(time));
    y_vec = zeros(1,length(time));
    
    
    for t = 1:length(time)
        
        x_dot = spike(t)*T_x_1*(1-x) - T_x_2*x;% +alpha*x_temp;
        x = x_dot/Fs + x;
        y_dot = (1-y)*x*T_y_1-y*T_y_2 + alpha*y; %*y^n/(y^n+k^n);
        y = y_dot/Fs + y;
        z = y^n/(y^n+k^n);
        
        x_vec(t) = x;
        y_vec(t) = y;
        z_vec(t) = z;
        
    end
    
    mean_exc(f) = mean(z_vec(2*Fs:3*Fs));
    p2p_exc(f) = max(z_vec(2*Fs:3*Fs))-min(z_vec(2*Fs:3*Fs));
    
    figure(1)
    plot(time,z_vec)
    hold on
    
end

fusion = 1-p2p_exc/p2p_exc(1);

figure(2)
plot(FR_test,mean_exc,'LineWidth',1)
xlabel('Frequency (Hz)','FontSize',14)
ylabel('Activation','FontSize',14)
hold on

figure(3)
plot(FR_test,fusion,'LineWidth',1)
xlabel('Frequency (Hz)','FontSize',14)
ylabel('Fusion','FontSize',14)
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
