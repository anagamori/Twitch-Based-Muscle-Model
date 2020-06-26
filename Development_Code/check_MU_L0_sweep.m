%==========================================================================
% check_MU_LO_sweep.m
% Author: Akira Nagamori
% Last update: 6/23/2020
% Descriptions: 
%   Check model characteristics (contraction time, twitch-tetanus ratio,
%   etc) of temporary model parameter sets generated by build_MU_LO_sweep
%==========================================================================

close all
clc
clear all

code_folder = '/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Development_Code';
data_folder = '/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Model Parameters/Model_11';

f_exp = [2 4 6 8 10 12.5 14 16 18 20 22 25 28 30 33.3 40 50 66.6 80 100];
force_exp = [2.5723 7.717 16.399 29.582 44.373 50.482 62.379 63.344 67.203 73.234 77.492 83.923 87.46 88.103 91.318 96.141 97.428 100.32 100 98.392];
[~,loc] = min(abs(force_exp-50));
f_half_exp = f_exp(loc);
fusion_exp = [-0.28169 17.465 36.62 58.31 74.93 85.915 90.141 92.394 94.93 96.056 96.62 97.465 97.465 97.746 98.028 97.746 97.465 97.465 97.746 96.62];


for i = 221
    i 
    cd('/Users/akiranagamori/Documents/Github/Twitch-Based-Muscle-Model/Development_Code/Data')
    load(['Data_temp_' num2str(i)])
    cd(code_folder)

    CT = Data{2,1}
    t2t = Data{2,5}
    
    figure()
    plot(Data{2,9},Data{2,11}*100)
    hold on
    plot(f_exp./f_half_exp,fusion_exp,'LineWidth',1,'color','k')
    xlim([0 3])
    %hold on 
    
end