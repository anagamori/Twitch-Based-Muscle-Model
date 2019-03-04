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

N_MU = 300; % number of motor units in a pool
load('index_slow') % index for the largest slow-twitch MU

parameter_Matrix = zeros(N_MU,15);
for i = 1:N_MU
    if i <= index_slow
        data_folder = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Data/ST';
    else
        data_folder = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Data/FT';
    end
    
    cd(data_folder)
    load(['MU_' num2str(i)])
    cd(code_folder)
    
    parameter_Matrix(i,:) = parameter;
    
end

save('pool_parameter_matrix','parameter_Matrix')