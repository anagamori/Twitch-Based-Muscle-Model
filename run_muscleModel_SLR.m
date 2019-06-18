%==========================================================================
% run_muscleModel.m
% Author: Akira Nagamori
% Last update: 3/5/19
% Descriptions:
%   Run muscle model simulation
%==========================================================================
close all
clear all
clc


%%
data_folder = '/Volumes/DATA2/New_Model/SLR/GD_40_GS_40_Ia_2000_DL_30';
code_folder = '/Users/akira/Documents/Github/Twitch-Based-Muscle-Model';
model_parameter_folder =  '/Users/akira/Documents/Github/Twitch-Based-Muscle-Model/Model Parameters/Model_4_Ur_50';

%%
cd(model_parameter_folder)
load('modelParameter')
cd(code_folder)
%% MU simulation parameters
modelParameter.CV_MU = 0.1;

%% Recruitment Type
modelParameter.recruitment = 3; % 1: Loeb's formulation, 2: Fuglevand's formulation


%% Simlulation parameters

amp_vec = [0.05 0.1:0.1:1];
trial_vec = [7 10];
for j = 1
    j
    if j < 2
        Fs = 10000;
        time = 0:1/Fs:15;
    elseif j >= 2 && j < 4
        Fs = 15000;
        time = 0:1/Fs:15;
    elseif j >= 4 && j < 7
        Fs = 20000;
        time = 0:1/Fs:15;
    elseif j >= 7
        Fs = 25000;
        time = 0:1/Fs:15;
    end
    %%
    SLRParameter.gamma_dynamic = 20;
    SLRParameter.gamma_static = 20;
    SLRParameter.Ia_delay = 30*Fs/1000;
    SLRParameter.Ia_gain = 2000;
    SLRParameter.Ib_gain = 2000;
    SLRParameter.Ib_delay = 40*Fs/1000;
    SLRParameter.RI_gain = 2;
    SLRParameter.RI_delay = 5*Fs/1000;

    amp = amp_vec(j+1);
    input = [zeros(1,1*Fs) amp/2*[0:1/Fs:2] amp*ones(1,length(time)-1*Fs-length(amp*[0:1/Fs:2]))];
    input_e =0.01*sin(2*pi*10*time);
    input(5*Fs+1:end) = input(5*Fs+1:end) + input_e(5*Fs+1:end);
    %%
    
    for i = 2
        i
        tic
        output = spikeDrivenMuscleModel_SLR(Fs,time,input,modelParameter,SLRParameter,1);
        toc
        cd(data_folder)
        save(['Data_' num2str(j) '_' num2str(i)],'output','-v7.3')
        cd(code_folder)
        clear output

    end
    
end