%==========================================================================
% test_MNs_v2.m
% Author: Akira Nagamori
% Last update: 6/21/19
%==========================================================================

close all
clear all
clc

code_folder = '/Users/akira/Documents/Github/Twitch-Based-Muscle-Model/Izhikevich Model';
model_parameter_folder =  '/Users/akira/Documents/Github/Twitch-Based-Muscle-Model/Model Parameters/Model_4_Ur_50_constantT2T';

%%
Fs = 10000;
time = 0:1/Fs:5;

%%
for n = 1:6
    cd(model_parameter_folder )
    load('modelParameter')
    cd(code_folder)
    FR_half = modelParameter.FR_half;
    MDR = modelParameter.FR_half/2;
    PDR = modelParameter.FR_half*2;
    U_th = modelParameter.U_th_new;
    testUnit = n;
    
    U_th_target = U_th(testUnit);
    MDR_target = MDR(testUnit);
    PDR_target = PDR(testUnit);
    FR_half_target = FR_half(testUnit);
    
    cd([model_parameter_folder '/MN'])
    load(['MN_' num2str(n)])
    cd(code_folder)
    
    input_vec = parameter.I_th-1:0.01:parameter.I_max;
    mean_FR = zeros(1,length(input_vec));
    for i = 1:length(input_vec)
        mean_FR(i) = FR_test(Fs,time,input_vec(i),parameter);
    end
 
    figure(1)
    subplot(2,1,1)
    plot(input_vec,mean_FR,'LineWidth',2)
    xlabel('Current (AU)','FontSize',14)
    ylabel('Discharge Rate (Hz)','FontSize',14)
    hold on
    plot([parameter.I_th-2 parameter.I_max+2],[MDR(testUnit) MDR(testUnit)],'k')
    plot([parameter.I_th-2 parameter.I_max+2],[PDR(testUnit) PDR(testUnit)],'k')
    plot([parameter.I_th parameter.I_th],[MDR(testUnit) PDR(testUnit)],'k')
    
    U_vec = 0:0.001:1;
    mean_FR = zeros(1,length(U_vec));
    for i = 1:length(U_vec)
        I = (parameter.I_max-parameter.I_th)/(1-U_th(testUnit))*(U_vec(i)-U_th(testUnit)) + parameter.I_th;
        mean_FR(i) = FR_test(Fs,time,I,parameter);
    end
    
    figure(1)
    subplot(2,1,2)
    plot(U_vec,mean_FR/FR_half_target,'LineWidth',2)
    xlabel('Synaptic Input (AU)','FontSize',14)
    ylabel('Discharge Rate (Hz)','FontSize',14)
    hold on
end