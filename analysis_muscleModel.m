%==========================================================================
% analysis_muscleModel.m
% Author: Akira Nagamori
% Last update: 3/5/19
% Descriptions:
%   Run muscle model simulation
%==========================================================================
close all
clear all
clc

%%
data_folder = '/Volumes/DATA2/New_Model/noTendon';
code_folder = '/Users/akira/Documents/Github/Twitch-Based-Muscle-Model';

%% 
amp_vec = 0.1:0.1:1;
Fs = 1000;

mean_Force = zeros(10,length(amp_vec));
std_Force = zeros(10,length(amp_vec));
cov_Force = zeros(10,length(amp_vec));
%% 
for j = 1:length(amp_vec)
    for i = 1:10
        load(['Data_' num2str(j) '_' num2str(trialN)])
        Force =  output.Force;
        mean_Force(i,j) = mean(Force(5*Fs+1:end));
        std_Force(i,j) = std(Force(5*Fs+1:end));
        cov_Force(i,j) =  std_Force(i,j)/meanForce(i,j);
    end
end

figure(1)
plot([0 amp_vec],[0 mean(mean_Force)])