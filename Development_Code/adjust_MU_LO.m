%==========================================================================
% adjust_MU_LO.m
% Author: Akira Nagamori
% Last update: 6/26/2020
% Descriptions:
%   Randomly adjust model parameters to generate temporary sets of
%   parameters
%==========================================================================


close all
clc
clear all

code_folder = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Development_Code';
data_folder = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Model Parameters/Model_11';

MU_type = 'fast';
Lce = 1;
FR_test_vec = [2:2:100 200 300];

MU_No = 1;
cd('/Users/akira/Documents/Github/Twitch-Based-Muscle-Model/Development_Code/Data')
load(['Data_v2_' num2str(MU_No)])
cd(code_folder)

CT = Data{2,1}
t2t = Data{2,5}
FR_half = Data{2,6}
    
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
phi_1 = parameter(13);
phi_2 = parameter(14);

a = -1;
b = 1;

r = a + (b-a)*rand(6,1); %0.3
r = r*0.1;
r_2 = a + (b-a)*rand(1,1); %0.4
r_2 = r_2*0.1;
%r_3 = rand(1,1)*0.1;

r_3 = a + (b-a)*rand(1,1);
r_3 = r_3 * 0; %0.1

S = S + S*r(1);
k_1 = k_1 + k_1*r(2);
k_2 = k_2 + k_2*r(3);
k_3 = k_3 + k_3*r(4);
k_4 = k_4 + k_4*r(5);
tau_1 = tau_1 - tau_1*r(6); %parameter(9); %0.0051;
K = K + K*r_2;
tau_3 = tau_3 + 0.03; %tau_3*r_3;


param = [S,C,k_1,k_2,k_3,k_4,tau_1,tau_2,N,K,tau_3,alpha,phi_1,phi_2];

simOpt = 0;
Fs = 2000;

[Data_new] = MUModel_test_full(param,Lce,0,FR_test_vec,MU_type,1,simOpt,Fs);

CT = Data_new{2,1}
t2t = Data_new{2,5}
FR_half = Data_new{2,6}
    
%%
Fs = 10000;
[Data_new] = MUModel_test_full(param,Lce,0,FR_test_vec,MU_type,1,simOpt,Fs);

%%
% cd('/Users/akira/Documents/Github/Twitch-Based-Muscle-Model/Development_Code/Data')
% save(['Data_v2_new_' num2str(MU_No)],'Data_new')
% cd(code_folder)
% clear Data

%%
Data = Data_new;
cd('/Users/akira/Documents/Github/Twitch-Based-Muscle-Model/Development_Code/Data')
save(['Data_v2_' num2str(MU_No)],'Data')
cd(code_folder)
