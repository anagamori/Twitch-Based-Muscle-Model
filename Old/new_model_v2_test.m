%==========================================================================
% new_model_v2_test.m
% Author: Akira Nagamori
% Last update: 11/9/18
%==========================================================================

close all
clear all
clc

%==========================================================================
% Simulation parameters
FR = 200;
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
x = 0; % proportion of Ca2+ bound to troponin (0-1)
y = 0; % porportion of cross-bridge sites available for binding (0-1)
z = 0; % porportion of bidning sites that formed cross-bridges (0-1)

T_x_1 = 400*Fs/1000; % forward rate constant for Ca2+ binding to troponin
T_x_2 = 50*Fs/1000; % backward rate constant for Ca2+ release from troponin 
T_y_1 = 200*Fs/1000; % forward rate constant for binding sites becoming available 
T_y_2 = 80*Fs/1000; % backward rate constant for binding sites becoming unavailable 
T_z_1 = 10*Fs/1000; % forward rate constant for binding sites forming cross-bridge 
T_z_2 = 40*Fs/1000; % backward rate constant for binding sites detachecing 

alpha = 10;

n = 2.5;
k = 0.3;
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
    y_dot = (1-y)*x*T_y_1-y*T_y_2 +alpha*z; %^n/(z^n+k^n);
    y = y_dot/Fs + y;
    y_int = y^n/(y^n+k^n);
    z_dot = (y_int-z)*T_z_1-z*T_z_2;
    z = z_dot/Fs + z;
    act = z; %z^n/(z^n+k^n);
    %R_vec(t) = R; 
    x_vec(t) = x;
    y_vec(t) = y;
    z_vec(t) = z;
    act_vec(t) = act;
end

if simulation_condition == 1
    % Twitch analysis
    [pks,locs_peak] = max(act_vec);
    pks
    CT = (locs_peak-1*Fs)*1000/Fs
    
    peak_half = pks/2;
    [~,HRT] = min(abs(act_vec(locs_peak:end)-peak_half));
    locs_hrt = locs_peak+HRT;
    HRT = HRT*1000/Fs
else
    mean_exc = mean(act_vec(2*Fs:3*Fs))
end

figure(1)
plot(time,x_vec)
hold on 
plot(time,y_vec)
legend('x','y')

figure(2)
plot(time,act_vec)
hold on
% plot(time,y_int_vec)

if simulation_condition == 1
    hold on
    plot([time(locs_peak) time(locs_peak)],[0 1],'--')
    text(time(locs_peak),pks,num2str(CT))
    hold on
    plot([time(locs_hrt) time(locs_hrt)],[0 1],'--')
    text(time(locs_hrt),peak_half,num2str(HRT))
    text(time(locs_peak)-time(locs_peak)*0.4,pks,num2str(pks))
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
