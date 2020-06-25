%==========================================================================
% build_MU_length.m
% Author: Akira Nagamori
% Last update: 6/23/2020
% Descriptions:
%   Randomly adjust model parameters to generate temporary sets of
%   parameters
%==========================================================================


close all
clc
clear all

code_folder = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Development_Code';
data_folder = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Model Parameters/Model_11';

MU_type = 'slow';

simOpt = 0;
Fs = 2000;


MU_No = 1;
cd('/Users/akira/Documents/Github/Twitch-Based-Muscle-Model/Development_Code/Data')
load(['Data_' num2str(MU_No)])
cd(code_folder)

parameter = Data{2,12};
S = parameter(1);
C = 1.8; %parameter(2);
k_1 = parameter(3);
k_2 = parameter(4);
k_3 = parameter(5);
k_4 = parameter(6);
tau_1 = parameter(7);
tau_2 = parameter(8);
N = parameter(9);
K = parameter(10);
tau_3 = parameter(11);
alpha = parameter(12);


Lce = 1;
[Data] = MUModel_test_L0(parameter,Lce,0,MU_type,1,simOpt,Fs);

hold on

FR_half = Data{2,6};
%%
a = 0;
b = 1;
r = a + (b-a)*rand(4,1); %0.3
r = r*0.3;
ptb = 0.1;

% k_3 = k_3 - k_3*r(1);
% k_4 = k_4 - k_4*r(2);
% N = N - N*r(3);
% K = K - K*r(4);

k_3 = k_3 - k_3*ptb; %r(1);
% k_4 = k_4 - k_4*ptb; %r(2);
N = N + N*ptb; %r(3);
alpha = alpha - alpha*ptb;
K = K + K*ptb; %r(4);


param = [S,C,k_1,k_2,k_3,k_4,tau_1,tau_2,N,K,tau_3,alpha];

Lce = 0.9;
[Data] = MUModel_test_L0(param,Lce,FR_half,MU_type,1,simOpt,Fs);

% cd('/Users/akira/Documents/Github/Twitch-Based-Muscle-Model/Development_Code/Data')
% save(['Data_temp_' num2str(MU_No)],'Data')
% cd(code_folder)
% clear Data

