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
%% Muscle architectural parameters
modelParameter.pennationAngle = 9.6*pi/180; %[radians]
modelParameter.optimalLength = 6.8; % [cm]
modelParameter.tendonSlackLength = 24.1; % [cm]
modelParameter.mass = 0.01; % [g]
modelParameter.muscleInitialLength = 6.8; % [cm]
modelParameter.tendonInitialLength = 24.1; % [cm]

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

%% Model parameters for activation-frequency relationship
load('pool_parameter_matrix')
modelParameter.parameterMatrix = parameter_Matrix;

%% FR_half for individual motor units
load('FR_half')
modelParameter.FR_half = FR_half;

%% Simlulation parameters
Fs = 5000;
time = 0:1/Fs:5;
amp = 0.1;
input = [zeros(1,1*Fs) amp/2*[0:1/Fs:2] amp*ones(1,length(time)-1*Fs-length(amp*[0:1/Fs:2]))];
%%
[output] = muscleModel_withTendon(Fs,time,input,modelParameter);