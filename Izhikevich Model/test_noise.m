%==========================================================================
% test_noise.m
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

[b,a] = butter(4,100/(Fs/2),'low');
%%
for n = 1
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
    
    U = 0.2; %0:0.001:1;
    r = normrnd(0,3,[1,length(time)]);
    r = filtfilt(b,a,r);
    I_amp = (parameter.I_max-parameter.I_th)/(1-U_th(testUnit))*(U-U_th(testUnit)) + parameter.I_th;
    input = [zeros(1,1*Fs) I_amp/2*[0:1/Fs:2] I_amp*ones(1,length(time)-1*Fs-length(I_amp*[0:1/Fs:2]))];
    input = input + r;
    [v_vec,binary] = Izhikevich(time,input,parameter,Fs);
    
    spike_time = find(binary(3*Fs+1:end));
    ISI = diff(spike_time)/(Fs/1000);
    FR = mean(1./ISI*1000)
    FR_sd = std(1./ISI*1000);
    ISI_CoV = FR_sd/FR*100
    
    figure(1)
    ax1 = subplot(2,1,1);
    plot(time,input,'LineWidth',2)
    ax2 = subplot(2,1,2);
    plot(time,v_vec,'LineWidth',2)
    linkaxes([ax1 ax2],'x')

end
