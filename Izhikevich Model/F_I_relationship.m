%==========================================================================
% F_I_relationship.m
% Author: Akira Nagamori
% Last update: 6/19/19
% Descriptions:
%   Test the frequency-current relationship of Izhikevich model
%==========================================================================

%close all
clear all
clc

code_folder = '/Users/akira/Documents/Github/Twitch-Based-Muscle-Model/Izhikevich Model';
model_parameter_folder =  '/Users/akira/Documents/Github/Twitch-Based-Muscle-Model/Model Parameters/Model_4_Ur_50_constantT2T';

%%
cd(model_parameter_folder)
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

input_vec = [0.1:0.1:0.9 1:100];
for i = 1:length(input_vec)
%% Input 
amp = input_vec(i);
input = [zeros(1,1*Fs) amp/2*[0:1/Fs:2] amp*ones(1,length(time)-1*Fs-length(amp*[0:1/Fs:2]))];

%% Model parameters
parameter.a = 0.04; %0.02 %the time scale of the recovery variable. Smaller values result in slower recovery. %0.005 for 2.4573 Hz 
parameter.b = 0.245; %0.2 %the sensitivity of the recovery variableu to the subthreshold fluctuations of the membrane potential v. 
%Greater values couple v and u more strongly resulting in possible subthreshold oscillations and low-threshold spiking dynamics
parameter.c = -65; % -65 %the after-spike reset value of the membrane potential v caused by the fast high-threshold K+ conductances
parameter.d = 200; %8 %after-spike reset of the recovery variable u caused by slowhigh-threshold Na+ andK+ conductances
parameter.v = -65; %-65

parameter.alpha = 0.04;
parameter.beta = 5;
parameter.gamma = 140;

%% Run Izhikevich model
[v_vec,binary] = Izhikevich(time,input,parameter,Fs);

spike_time = find(binary(3*Fs+1:end));
ISI = diff(spike_time)/(Fs/1000);
mean_FR(i) = mean(1./ISI*1000);

end

figure(1)
plot(input_vec,mean_FR)
hold on
plot([0 100],[MDR(testUnit) MDR(testUnit)],'k')
plot([0 100],[PDR(testUnit) PDR(testUnit)],'k')