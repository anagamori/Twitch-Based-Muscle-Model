close all
clc
clear all

code_folder = '/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Development_Code';
data_folder = '/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Model Parameters/Model_11';

MU_type = 'fast';
Lce = 1;
S = 17.3; %parameter(1); %7;
C = 1.8; %parameter(2); %1.025;
k_1 = 76; %parameter(3); %14.625;
k_2 = 1.1; %parameter(4); %4.9375;
k_3 = 59; %parameter(5)*Lce + parameter(6); %17.41*Lce - 2.85;
k_4 = 38; %parameter(7)*Lce + parameter(8); %-7.67*Lce + 14.92;
tau_1 = 0.002; %parameter(9); %0.0051;
tau_2 = 0.025; %parameter(10); % 0.04;
N = 2.2; %parameter(11)*Lce + parameter(12); %-2.26*Lce + 4.20;
K = 0.0623; %parameter(13)*Lce + parameter(14); %-0.044*Lce + 0.080;
tau_3 = 0.1; %parameter(15); %4.475;
alpha = 0.8;

param = [S,C,k_1,k_2,k_3,k_4,tau_1,tau_2,N,K,tau_3,alpha];

simOpt = 1;
Fs = 2000;

[~] = MUModel_test_v3(param,Lce,0,MU_type,1,simOpt,Fs);
