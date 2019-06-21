%==========================================================================
% build_MNs.m
% Author: Akira Nagamori
% Last update: 6/18/19
%==========================================================================


close all
clear all
clc

code_folder = '/Users/akiranagamori/Documents/Github/Twitch-Based-Muscle-Model/Izhikevich Model';
model_parameter_folder =  '/Users/akiranagamori/Documents/Github/Twitch-Based-Muscle-Model/Model Parameters/Model_4_Ur_50';

%%
cd(model_parameter_folder )
load('modelParameter')
cd(code_folder)
MDR = modelParameter.FR_half/2;
PDR = modelParameter.FR_half*2;
U_th = modelParameter.U_th_new;
[val,loc] = min(U_th);
testUnit = loc;
%%
Fs = 10000;
time = 0:1/Fs:5;

%% Input 
%input = zeros(1,length(time));
% input = zeros(1,length(time));
% input(2*Fs+1:end) = 1;
amp = 10;
input = [zeros(1,1*Fs) amp/2*[0:1/Fs:2] amp*ones(1,length(time)-1*Fs-length(amp*[0:1/Fs:2]))];
%% Model parameters
parameter.a = 0.02; %; %0.02 %the time scale of the recovery variable. Smaller values result in slower recovery. %0.005 for 2.4573 Hz 
parameter.b = 0.2; %0.2 %the sensitivity of the recovery variableu to the subthreshold fluctuations of the membrane potential v. 
%Greater values couple v and u more strongly resulting in possible subthreshold oscillations and low-threshold spiking dynamics
parameter.c = -65; % -65 %the after-spike reset value of the membrane potential v caused by the fast high-threshold K+ conductances
parameter.d = 10; %8 %after-spike reset of the recovery variable u caused by slowhigh-threshold Na+ andK+ conductances
parameter.v = -64; %-65

parameter.alpha = 0.04;
parameter.beta = 5;
parameter.gamma = 140;

%% Run Izhikevich model
[v_vec,binary] = Izhikevich(time,input,parameter,Fs);

figure(1)
ax1 = subplot(2,1,1);
plot(time,input)
ax2 = subplot(2,1,2);
plot(time,v_vec)
linkaxes([ax1 ax2],'x')

spike_time = find(binary(3*Fs+1:end));
ISI = diff(spike_time)/(Fs/1000);
mean_FR = mean(1./ISI*1000)
minimum_rate = MDR(testUnit)
maximum_rate = PDR(testUnit)