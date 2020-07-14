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
data_folder = '/Volumes/DATA2/PLOS_CB_Data/withTendon/onion_skin';
code_folder = '/Users/akira/Documents/Github/Twitch-Based-Muscle-Model';
model_parameter_folder =  '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Development_Code/Data';

%%
cd(model_parameter_folder)
load('modelParameter_onion_skin')
cd(code_folder)

%% Simlulation parameters

amp_vec = [0.025 0.05 0.1:0.1:1];
%amp_vec = [0.106 0.29 0.62 0.78 0.93];
trial_vec = [7 10];
for j = 1:length(amp_vec)
    j
%      Fs = 20000;
%         time = 0:1/Fs:15;
    if j <= 7
        Fs = 10000;
        time = 0:1/Fs:15;
    elseif j > 7 && j <= 9
        Fs = 15000;
        time = 0:1/Fs:15;
    elseif j >= 10
        Fs = 20000;
        time = 0:1/Fs:15;
    end
   
    amp = amp_vec(j);
    input = [zeros(1,1*Fs) amp/2*[0:1/Fs:2] amp*ones(1,length(time)-1*Fs-length(amp*[0:1/Fs:2]))];
    %%
    for i = 7:10
        i
        tic
        output = MU_population_model_onion_skin(Fs,time,input,modelParameter,1);
        toc
        cd(data_folder)
        save(['Data_' num2str(j) '_' num2str(i)],'output','-v7.3')
        cd(code_folder)
        Force = output.ForceTendon(end-5*Fs+1:end);
        mean(Force(end-1*Fs+1:end))/444.5844
        CoV = std(Force)/mean(Force)
        clear output
    end

end

