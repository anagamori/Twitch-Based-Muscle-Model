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
data_folder = '/Volumes/DATA2/PLOS_CB_Data/noTendon/Model_onion_skin_N_100_Ur_50';
code_folder = '/Users/akira/Documents/Github/Twitch-Based-Muscle-Model';
model_parameter_folder =  '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Development_Code/Data';

%%
cd(model_parameter_folder)
load('modelParameter_onion_skin_N_100_Ur_50')
cd(code_folder)

%% Simlulation parameters

amp_vec = [0.025 0.05 0.1:0.1:1];
%amp_vec = 0:0.01:1;
trial_vec = [7 10];
for j = 6:length(amp_vec)
    j

        Fs = 10000;
        time = 0:1/Fs:15;

    amp = amp_vec(j);
    input = [zeros(1,1*Fs) amp/2*[0:1/Fs:2] amp*ones(1,length(time)-1*Fs-length(amp*[0:1/Fs:2]))];
 
    %%
    for i = 3:4
        i
        tic
        output = MU_population_model_noTendon(Fs,time,input,modelParameter,1);
        toc
        cd(data_folder)
        save(['Data_' num2str(j) '_' num2str(i)],'output','-v7.3')
        cd(code_folder)
        Force = output.Force(end-5*Fs+1:end);
        CoV = std(Force)/mean(Force)
        clear output
    end
 
    
    
end

