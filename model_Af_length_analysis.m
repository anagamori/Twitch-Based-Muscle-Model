%==========================================================================
% new_model_parameterFit_ST_length.m
% Author: Akira Nagamori
% Last update: 2/22/119
%==========================================================================

close all
clear all
clc

%% Folder name
code_folder = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model';
data_folder = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Data';

%% 
trialN = 1;
param_Matrix_1 = zeros(10,11);
param_Matrix_2 = zeros(10,11);
param_Matrix_3 = zeros(10,11);
param_Matrix_4 = zeros(10,11);
for i = 1:10
    
    cd(data_folder)
    load(['Data_' num2str(trialN) '_' num2str(i+10)])
    cd(code_folder)    
    param_Matrix_1(i,:) = Data{2,12};
    
    cd(data_folder)
    load(['Data_' num2str(trialN) '_' num2str(i+20)])
    cd(code_folder)
    param_Matrix_2(i,:) = Data{2,12};
    
    cd(data_folder)
    load(['Data_' num2str(trialN) '_' num2str(i+30)])
    cd(code_folder)
    param_Matrix_3(i,:) = Data{2,12};
    
    cd(data_folder)
    load(['Data_' num2str(trialN) '_' num2str(i+40)])
    cd(code_folder)
    param_Matrix_4(i,:) = Data{2,12};
    
end