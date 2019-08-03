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
data_folder = '/Volumes/DATA2/New_Model/withTendon/Model_6_20_CoV_50_Ur_Rec_3_noFV';
code_folder = '/Users/akira/Documents/Github/Twitch-Based-Muscle-Model';
model_parameter_folder =  '/Users/akira/Documents/Github/Twitch-Based-Muscle-Model/Model Parameters/Model_6';

%%
cd(model_parameter_folder)
load('modelParameter')
cd(code_folder)
%% MU simulation parameters
modelParameter.CV_MU = 0.2;

%% Recruitment Type
modelParameter.recruitment = 3; % 1: Loeb's formulation, 2: Fuglevand's formulation


%% Simlulation parameters

amp_vec = [0.05 0.1:0.1:1];
trial_vec = [7 10];
for j = 0:10
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
        elseif j >= 7 && j <= 8
            Fs = 25000;
            time = 0:1/Fs:15;
        elseif j >= 9
            Fs = 30000;
            time = 0:1/Fs:15;
        end
%     if j <= 1
%         Fs = 10000;
%         time = 0:1/Fs:15;
%     elseif j >= 2 && j < 3 
%         Fs = 15000;
%         time = 0:1/Fs:15;
%     elseif j >= 3 && j <= 4
%         Fs = 20000;
%         time = 0:1/Fs:15;
%     elseif j >= 5 && j < 6
%         Fs = 25000;
%         time = 0:1/Fs:15;
%     elseif j >= 6 && j <= 7 
%         Fs = 30000;
%         time = 0:1/Fs:15;
%     else
%         Fs = 40000;
%         time = 0:1/Fs:15;
%     end
    amp = amp_vec(j+1);
    input = [zeros(1,1*Fs) amp/2*[0:1/Fs:2] amp*ones(1,length(time)-1*Fs-length(amp*[0:1/Fs:2]))];
    %%
    if j == 10
        for i = 1:10
            i
            tic
            output = spikeDrivenMuscleModel(Fs,time,input,modelParameter,1);
            toc
            cd(data_folder)
            save(['Data_' num2str(j) '_' num2str(i)],'output','-v7.3')
            cd(code_folder)
            clear output
            
        end
    else
        for i = 1:10
            i
            tic
            output = spikeDrivenMuscleModel(Fs,time,input,modelParameter,1);
            toc
            cd(data_folder)
            save(['Data_' num2str(j) '_' num2str(i)],'output','-v7.3')
            cd(code_folder)
            clear output
            
        end
    end
end