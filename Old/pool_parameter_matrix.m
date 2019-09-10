%==========================================================================
% pool_parameter_matrix.m
% Author: Akira Nagamori
% Last update: 3/4/19
% Descriptions: Create a matrix of parameters that define the
% activation-frequency relationship of each motor unit
%==========================================================================

close all
clc
clear all

code_folder = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model';
data_folder_default = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Model Parameters/Model_2/';

N_MU = 300; % number of motor units in a pool
cd(data_folder_default)
load('index_slow') % index for the largest slow-twitch MU
cd(code_folder)
parameter_Matrix = zeros(N_MU,15);
for i = 1:N_MU
    if i <= index_slow
        data_folder = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Model Parameters/Model_2/ST';
    else
        data_folder = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Model Parameters/Model_2/FT';
    end
    
    cd(data_folder)
    load(['MU_' num2str(i)])
    cd(code_folder)
    
    parameter_Matrix(i,:) = parameter;
    
end

cd(data_folder_default)
save('pool_parameter_matrix','parameter_Matrix')
cd(code_folder)