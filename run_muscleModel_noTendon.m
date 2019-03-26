%==========================================================================
% run_muscleModel_noTendon.m
% Author: Akira Nagamori
% Last update: 3/5/19
% Descriptions:
%   Run muscle model simulation
%==========================================================================
close all
clear all
clc

%%
data_folder = '/Volumes/DATA2/New_Model/noTendon/No_100_adjustedForce';
code_folder = '/Users/akira/Documents/Github/Twitch-Based-Muscle-Model';
model_parameter_folder =  '/Users/akira/Documents/Github/Twitch-Based-Muscle-Model/Model Parameters/Model_N_100_adjustedF0';

%% 
cd(model_parameter_folder)
load('modelParameter')
cd(code_folder)
%% MU simulation parameters
modelParameter.CV_MU = 0.1;

%% Recruitment Type
modelParameter.recruitment = 2; % 1: Loeb's formulation, 2: Fuglevand's formulation


%% Simlulation parameters
Fs = 1000;
time = 0:1/Fs:15;
amp_vec = 0.1:0.1:1;
for j = 1:length(amp_vec)
    j
    amp = amp_vec(j);
    input = [zeros(1,1*Fs) amp/2*[0:1/Fs:2] amp*ones(1,length(time)-1*Fs-length(amp*[0:1/Fs:2]))];
    
    output_temp = cell(1,10);
    
    parfor i = 1:10
        tic
        output_temp{i} = muscleModel_noTendon(Fs,time,input,modelParameter);
        toc
        
    end
    for trialN = 1:10
        output = output_temp{trialN};
        cd(data_folder)
        save(['Data_' num2str(j) '_' num2str(trialN)],'output')
        cd(code_folder)
    end
%     for i = 1:10
%         tic
%         output = muscleModel_noTendon(Fs,time,input,modelParameter);
%         toc
%         cd(data_folder)
%         save(['Data_' num2str(j) '_' num2str(i)],'output')
%         cd(code_folder)
%     end
end