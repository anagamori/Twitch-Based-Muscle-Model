close all
clc
clear all
%%
code_folder = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Development_Code';
data_folder = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Model Parameters/Model_11';

%% 
% cd(data_folder)
% load('modelParameter')
% cd(code_folder)
%parameterMatrix = modelParameter.parameterMatrix;

cd(data_folder)
load('pool_parameter_matrix')
cd(code_folder)


% up to N = 147 => slow
MU_No = 100;
MU_type = 'slow';
Lce = 1;

%[S,C,k_1,k_2,k_3,k_4,tau_1,tau_2,N,K,tau_3]
parameter_temp = parameter_Matrix(MU_No,:);
parameter = parameter_temp;
parameter(15) = parameter(10);
parameter(10) = 0.02;

S = 18; %parameter(1); %7;
C = 1.8; %parameter(2); %1.025;
k_1 = 55; %parameter(3); %14.625;
k_2 = 50; %parameter(4); %4.9375;
k_3 = 30; %parameter(5)*Lce + parameter(6); %17.41*Lce - 2.85;
k_4 = 25; %parameter(7)*Lce + parameter(8); %-7.67*Lce + 14.92;
tau_1 = 0.003; %parameter(9); %0.0051;
tau_2 = 0.015; %parameter(10); % 0.04;
N = 2.2; %parameter(11)*Lce + parameter(12); %-2.26*Lce + 4.20;
K = 0.04; %parameter(13)*Lce + parameter(14); %-0.044*Lce + 0.080;
tau_3 = 0.04; %parameter(15); %4.475;
alpha = 1.2;

param = [S,C,k_1,k_2,k_3,k_4,tau_1,tau_2,N,K,tau_3,alpha];

simOpt = 0;
Fs = 2000;

[Data_temp] = MUModel_test_v2(param,Lce,0,MU_type,1,simOpt,Fs);
parameter = param;
%%
% cd('/Users/akira/Documents/Github/Twitch-Based-Muscle-Model/Development_Code/Data')
% save('seed_ST_1','parameter')
% cd(code_folder)