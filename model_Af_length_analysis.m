%==========================================================================
% new_model_parameterFit_ST_length.m
% Author: Akira Nagamori
% Last update: 2/22/119
%==========================================================================

close all
clear all
clc

%% Folder name
code_folder = '/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model';
data_folder = '/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Data';

%% 
trialN = 1;
param_Matrix_1 = zeros(10,11);
param_Matrix_2 = zeros(10,11);
param_Matrix_3 = zeros(10,11);
param_Matrix_4 = zeros(10,11);

cd(data_folder)
load(['Data_' num2str(trialN)])
param = Data{2,12};
cd(code_folder)
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

%% 
L0_long = [0.8 0.9 1.1 1.2];
L0_vec = [0.8*ones(10,1);0.9*ones(10,1);1.1*ones(10,1);1.2*ones(10,1)];

k_3_long = [param_Matrix_1(:,5) param_Matrix_2(:,5) param_Matrix_3(:,5) param_Matrix_4(:,5)]';
k_3_vec = reshape(k_3_long',[],1);

k_4_long = [param_Matrix_1(:,6) param_Matrix_2(:,6) param_Matrix_3(:,6) param_Matrix_4(:,6)]';
k_4_vec = reshape(k_4_long',[],1);

tau_1_long = [param_Matrix_1(:,7) param_Matrix_2(:,7) param_Matrix_3(:,7) param_Matrix_4(:,7)]';
tau_1_vec = reshape(tau_1_long',[],1);



