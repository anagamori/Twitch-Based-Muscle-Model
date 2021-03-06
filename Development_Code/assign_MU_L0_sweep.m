%==========================================================================
% assign_MU_LO_sweep.m
% Author: Akira Nagamori
% Last update: 6/23/2020
% Descriptions: 
%   Assign model parameters generated by build_MU_LO_sweep.m  to individual motor units based on their
%   contraction time 
%==========================================================================


close all
clc
clear all

code_folder = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Development_Code';
data_folder = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Model Parameters/Model_11';

% 222
for i = 233
    MU_No = 50;
    
    i 
    cd('/Users/akira/Documents/Github/Twitch-Based-Muscle-Model/Development_Code/Data')
    load(['Data_temp_' num2str(i)]) 
    cd(code_folder)

    CT = Data{2,1}
    t2t = Data{2,5}
    DR_half = Data{2,6}
    
    figure()
    plot(Data{2,9},Data{2,11})
    xlim([0 3])
    
    
    
    
    cd('/Users/akira/Documents/Github/Twitch-Based-Muscle-Model/Development_Code/Data')
    save(['Data_' num2str(MU_No)],'Data')
    cd(code_folder)

end