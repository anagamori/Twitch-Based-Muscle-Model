%==========================================================================
% new_model_v11_parameterSearch_visualization.m
% Author: Akira Nagamori
% Last update: 2/22/119
% Descriptions
%   Inspired by Williams et al. 1998
%==========================================================================

close all
clear all
clc

code_folder = '/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model';
data_folder = '/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Data';
    
%%
j = 1;
for m = 1501:2000
    cd (data_folder)
    load(['Data_' num2str(m)])   
    cd(code_folder)
    error_long(j) = Data{2,8};
    j = j + 1;
end

%%
index = find(error_long<0.4);
for i = 1:length(index)
    cd (data_folder)
    load(['Data_' num2str(index(i))])   
    cd(code_folder)
    
    CT_long(i) =  Data{2,1};
end

CT_long