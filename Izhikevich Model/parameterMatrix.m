%==========================================================================
% parameterMatrix.m
% Author: Akira Nagamori
% Last update: 6/21/19
%==========================================================================

close all
clear all
clc

code_folder = '/Users/akiranagamori/Documents/Github/Twitch-Based-Muscle-Model/Izhikevich Model';
model_parameter_folder =  '/Users/akiranagamori/Documents/Github/Twitch-Based-Muscle-Model/Model Parameters/Model_4_Ur_50_constantT2T';

%%
Fs = 10000;
time = 0:1/Fs:5;

a = zeros(1,300);
I_th = zeros(1,300);
I_max = zeros(1,300);

%%
for n = 1:300
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
    
    I_th(n) = parameter.I_th;
    I_max(n) = parameter.I_max;
    a(n) = parameter.a;
    
    
end

parameterMN.a = a;
parameterMN.I_th = I_th;
parameterMN.I_max = I_max;

cd(model_parameter_folder)
save('parameterMN','parameterMN')
cd(code_folder)
