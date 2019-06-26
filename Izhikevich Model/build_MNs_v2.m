%==========================================================================
% build_MNs_v2.m
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
for n = 1
    cd(model_parameter_folder )
    load('modelParameter')
    cd(code_folder)
    MDR = modelParameter.FR_half/2;
    PDR = modelParameter.FR_half*2;
    U_th = modelParameter.U_th_new;
    testUnit = n;
    
    U_th_target = U_th(testUnit);
    MDR_target = MDR(testUnit);
    PDR_target = PDR(testUnit);
    
    %% Fixed parameters
    parameter.a = 0.02;
    parameter.b = 0.2;
    parameter.c = -65;
    parameter.d = 6;
    parameter.v = -65;
    parameter.alpha = 0.04;
    parameter.beta = 5;
    parameter.gamma = 140;
    
    %% Find a threshold current, I_th
    input_vec = 0.5:0.01:50; %[0.5:0.1:0.9 1:100];
    mean_FR = zeros(1,length(input_vec));
    for i = 1:length(input_vec)
        mean_FR(i) = FR_test(Fs,time,input_vec(i),parameter);
    end
    index_nan = isnan(mean_FR);
    parameter.I_th = input_vec(find(index_nan==1,1,'last')+1);
    
    %% Find a parameter value of a that produces minimun discharge rate at I_th
    a_vec = 0:0.001:0.1;
    input_vec = parameter.I_th-1:0.01:parameter.I_th+1;
    FR_mat = zeros(length(a_vec),length(input_vec));
    for j = 1:length(a_vec)
        parameter.a = a_vec(j);
        for k = 1:length(input_vec)
        FR_mat(j,k) = FR_test(Fs,time,input_vec(k),parameter);
        end
    end
    FR_diff_mat = abs(MDR_target-FR_mat);
    FR_th = zeros(1,length(a_vec));
    index_th = zeros(1,length(a_vec));
    for k = 1:length(a_vec)
        temp = isnan(FR_mat(k,:));
        index_1 = find(temp==1,1,'last')+1;
        index_th(k) = index_1;
        if index_1 > length(input_vec)
            FR_th(k) = 0;
        else            
            FR_th(k) = FR_mat(k,index_1);
        end
    end
    [error,a_index] = min(abs(MDR_target-FR_th));
    parameter.a = a_vec(a_index);
    parameter.I_th = input_vec(index_th(a_index));
    
    %% Re-test frequency-current relationship and find I_min and I_max
    % I_max: current required to reach peak discharge rate
    input_vec = 0.5:0.01:50; %[0.5:0.1:0.9 1:100];
    mean_FR = zeros(1,length(input_vec));
    for i = 1:length(input_vec)
        mean_FR(i) = FR_test(Fs,time,input_vec(i),parameter);
    end
        
    [error_I_max,I_max_index] = min(abs(PDR_target-mean_FR));
    parameter.I_max = input_vec(I_max_index);
    
    %% Plot results
    figure(1)
    plot(input_vec,mean_FR)
    hold on
    plot([0 100],[MDR(testUnit) MDR(testUnit)],'k')
    plot([0 100],[PDR(testUnit) PDR(testUnit)],'k')
    plot([U_th_target U_th_target],[MDR(testUnit) PDR(testUnit)],'k')
    
    %% Save parameters
    cd([model_parameter_folder '/MN'])
    save(['MN_' num2str(n)],'parameter')
    cd(code_folder)
end