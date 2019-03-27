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
data_folder = '/Volumes/DATA2/New_Model/withTendon/10_CoV_50_Ur_Rec_2_shortTendon';
code_folder = '/Users/akira/Documents/Github/Twitch-Based-Muscle-Model';
%% Muscle architectural parameters
modelParameter.pennationAngle = 9.6*pi/180; %[radians]
modelParameter.optimalLength = 6.8; % [cm]
modelParameter.tendonSlackLength = 12.1; % [cm]
modelParameter.mass = 0.01; % [g]
modelParameter.muscleInitialLength = 6.8; % [cm]
modelParameter.tendonInitialLength = 12.1; % [cm]

%% MU simulation parameters
modelParameter.CV_MU = 0.1;
%% Contraction time
% Generate a distribution of contraction time across motor units based on
% Rayleigh distribution
% rng(1)
% min_CT = 20; %minimum contraction time [ms]
% CT = round(raylrnd(23,1,N_MU)+min_CT); %contraction time of individual motor units [ms]
% CT_sorted = sort(CT,'descend');
load('CT_vec')
modelParameter.CT = CT_vec;

%% Recruitment threshold
modelParameter.Ur = 0.5;

%% Range of peak tetanic tension
modelParameter.RP = 25;

%% Model parameters for activation-frequency relationship
load('pool_parameter_matrix')
modelParameter.parameterMatrix = parameter_Matrix;

%% FR_half for individual motor units
load('FR_half')
modelParameter.FR_half = FR_half;

%% Recruitment Type
modelParameter.recruitment = 2; % 1: Loeb's formulation, 2: Fuglevand's formulation

%% Simlulation parameters

amp_vec = 0.1:0.1:1;
trial_vec = [7 10];
for j = 7 %:length(amp_vec)
    j
    if j <= 2
        Fs = 10000;
        time = 0:1/Fs:15;
    elseif j > 2 && j <= 3
        Fs = 15000;
        time = 0:1/Fs:15;
    elseif j >= 4 && j <6
        Fs = 20000;
        time = 0:1/Fs:15;
    elseif j >= 6
        Fs = 25000;
        time = 0:1/Fs:15;
    end
    amp = amp_vec(j);
    input = [zeros(1,1*Fs) amp/2*[0:1/Fs:2] amp*ones(1,length(time)-1*Fs-length(amp*[0:1/Fs:2]))];
    %%
    
    
    for i = 1:4
        tic
        output = muscleModel_withTendon(Fs,time,input,modelParameter);
        toc
        cd(data_folder)
        save(['Data_' num2str(j) '_' num2str(i)],'output')
        cd(code_folder)
    end
    
end