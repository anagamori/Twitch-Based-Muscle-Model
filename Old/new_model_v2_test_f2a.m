%==========================================================================
% new_model_v2_test_f2a.m
% Author: Akira Nagamori
% Last update: 10/31/18
%
%==========================================================================

%close all
clear all
clc

%==========================================================================
% Simulation parameters
FR_test = [2 5 10 15 20 30 50 100 200]; %10:10:100];
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
    y = 0; % porportion of cross-bridge sites available for binding (0-1)
    z = 0; % porportion of bidning sites that formed cross-bridges (0-1)
    act = 0;
    T_x_1 = 400*Fs/1000; % forward rate constant for Ca2+ binding to troponin
    T_x_2 = 50*Fs/1000; % backward rate constant for Ca2+ release from troponin
    T_y_1 = 200*Fs/1000; % forward rate constant for binding sites becoming available
    T_y_2 = 80*Fs/1000; % backward rate constant for binding sites becoming unavailable
    T_z_1 = 10*Fs/1000; % forward rate constant for binding sites forming cross-bridge
    T_z_2 = 40*Fs/1000; % backward rate constant for binding sites detachecing

    alpha = 10;
    
    n = 3;
    k = 0.25;
    
    N = 2;
    %==========================================================================
    % initialization
    R = zeros(1,length(time));
    x_1_int_vec = zeros(1,length(time));
    x_vec = zeros(1,length(time));
    y_vec = zeros(1,length(time));
    y_int_vec = zeros(1,length(time));
    z_vec = zeros(1,length(time));
    act_vec = zeros(1,length(time));
    
    for t = 1:length(time)
        x_dot = spike(t)*(1-x)*T_x_1 - T_x_2*x;
        x = x_dot/Fs + x;
        y_dot = (1-y)*x*T_y_1-y*T_y_2 + alpha*act^3; %^n/(z^n+k^n);
        y = y_dot/Fs + y;
        y_int = y^n/(y^n+k^n);
        z_dot = (y_int-z)*T_z_1-z*T_z_2;
        z = z_dot/Fs + z;
        act = z*(T_z_1+T_z_2)/(T_z_1); %z^n/(z^n+k^n);
        %R_vec(t) = R;
        x_vec(t) = x;
        y_vec(t) = y;
        y_int_vec(t) = y_int;
        z_vec(t) = z;
        act_vec(t) = act;
        
    end
    
    mean_exc(f) = mean(act_vec(2*Fs:3*Fs));
    p2p_exc(f) = max(act_vec(2*Fs:3*Fs))-min(act_vec(2*Fs:3*Fs));
    
    
    figure(1)
    plot(time,act_vec)
    hold on 
%     hold on
%     plot(time,x_vec)
%     hold on 
%     plot(time,y_vec)
%     hold on 
%     plot(time,y_int_vec)
    
end

twitch2tetanus_ratio = p2p_exc(1)/mean_exc(f)
fusion = 1-p2p_exc/p2p_exc(1);

f_env = 0:0.01:FR_test(end)/22;
Af = 1-exp(-(f_env/(0.56*2.11)).^2.11);

figure(21)
plot(FR_test/17.5,mean_exc,'LineWidth',1)
hold on
plot(f_env,Af)
xlabel('Frequency (Hz)','FontSize',14)
ylabel('Activation','FontSize',14)
hold on
xlim([0 3])

figure(22)
plot(FR_test,fusion,'LineWidth',1)
xlabel('Frequency (Hz)','FontSize',14)
ylabel('Fusion','FontSize',14)
hold on

figure(23)
plot(mean_exc./max(mean_exc),fusion,'LineWidth',1)
xlabel('Activation','FontSize',14)
ylabel('Fusion','FontSize',14)
hold on
plot(0:0.1:1,0:0.1:1,'--')

function spikeTrain = spikeTrainGenerator(t,Fs,freq)

spikeTrain = zeros(1,length(t));
ISI = round(1/freq*Fs);
numSpikes = round(length(t)/ISI);
index = [1:numSpikes]*ISI;
index(index>length(t)) = [];
spikeTrain(index) = 1;
spikeTrain(1) = 1;

end
