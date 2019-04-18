%==========================================================================
% pool_FR_half.m
% Author: Akira Nagamori
% Last update: 3/4/19
% Descriptions: Create a vector of FR_half by running the model with the
% parameter set (pool_parameter_matrix)
%==========================================================================

close all
clc
clear all

code_folder = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model';
data_folder_default = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Model Parameters/Model_2/';

N_MU = 10; % number of motor units in a pool
cd(data_folder_default)
load('index_slow') % index for the largest slow-twitch MU
cd(code_folder)
parameter_Matrix = zeros(N_MU,15);
for i = 1:N_MU
    if i <= index_slow
        data_folder = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Model Parameters/Model_2/ST';
        MU_type = 'slow';
    else
        data_folder = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Model Parameters/Model_2/FT';
        MU_type = 'fast';
    end
    
    cd(data_folder)
    load(['MU_' num2str(i)])
    cd(code_folder)
    
    [Data] = spikeDrivenMuscleModel_testFunction_fullVersion(parameter,1,0,MU_type,1);
    FR_half(i) = Data{2,6};
    t2t(i) = Data{2,5};
    CT_vec(i) = Data{2,1};
    twitch_force(i) = Data{2,4};
    
end

%% 
cd(data_folder_default)
save('CT_vec','CT_vec')
save('FR_half','FR_half')
save('t2t','t2t')
save('twitch_force','twitch_force')
cd(code_folder)
MFR_MU = FR_half./2;
PFR_MU = FR_half.*2;