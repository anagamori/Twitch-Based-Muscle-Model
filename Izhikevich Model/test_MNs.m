%==========================================================================
% test_MNs.m
% Author: Akira Nagamori
% Last update: 6/18/19
%==========================================================================

close all
clear all
clc

code_folder = '/Users/akiranagamori/Documents/Github/Twitch-Based-Muscle-Model/Izhikevich Model';
model_parameter_folder =  '/Users/akiranagamori/Documents/Github/Twitch-Based-Muscle-Model/Model Parameters/Model_4_Ur_50_constantT2T';

%%
Fs = 10000;
time = 0:1/Fs:5;

%%
for n = 86
    cd(model_parameter_folder )
    load('modelParameter')
    cd(code_folder)
    MDR = modelParameter.FR_half/2;
    PDR = modelParameter.FR_half*2;
    U_th = modelParameter.U_th_new;
    testUnit = n;
    
    U_th_target = U_th(testUnit)*100;
    MDR_target = MDR(testUnit);
    PDR_target = PDR(testUnit);
    
    cd([model_parameter_folder '/MN'])
    load(['MN_' num2str(n)])
    cd(code_folder)
    
    input_vec = parameter.I_min-1:0.1:parameter.I_max;
     mean_FR = zeros(1,length(input_vec));
    for i = 1:length(input_vec)
        %% Input
        input = zeros(1,length(time));
        input(2*Fs+1:end) = input_vec(i);
        %% Run Izhikevich model
        [v_vec,binary] = Izhikevich(time,input,parameter,Fs);
        
        spike_time = find(binary(3*Fs+1:end));
        ISI = diff(spike_time)/(Fs/1000);
        mean_FR(i) = mean(1./ISI*1000);
        
    end
 
    figure(1)
    subplot(2,1,1)
    plot(input_vec,mean_FR,'LineWidth',2)
    xlim([parameter.I_min-1 parameter.I_max+1])
    xlabel('Current (AU)','FontSize',14)
    ylabel('Discharge Rate (Hz)','FontSize',14)
    
    U_vec = 0:0.001:1;
    mean_FR = zeros(1,length(U_vec));
    for i = 1:length(U_vec)
        I = (parameter.I_max-parameter.I_min)/(1-U_th(testUnit))*(U_vec(i)-U_th(testUnit)) + parameter.I_min;
        input = zeros(1,length(time));
        input(2*Fs+1:end) = I;
        
        %% Run Izhikevich model
        [v_vec,binary] = Izhikevich(time,input,parameter,Fs);
        
        spike_time = find(binary(3*Fs+1:end));
        ISI = diff(spike_time)/(Fs/1000);
        mean_FR(i) = mean(1./ISI*1000);
    end
    
    figure(1)
    subplot(2,1,2)
    plot(U_vec,mean_FR,'LineWidth',2)
    xlabel('Synaptic Input (AU)','FontSize',14)
    ylabel('Discharge Rate (Hz)','FontSize',14)
end