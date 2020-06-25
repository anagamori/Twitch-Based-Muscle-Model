%==========================================================================
% build_MU_LO_sweep.m
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

MU_type = 'fast';
Lce = 1;



%
% S = 14.10566883; %parameter(1); %7;
% C = 1.8; %parameter(2); %1.025;
% k_1 = 52.96651207; %parameter(3); %14.625;
% k_2 = 42.57566516; %parameter(4); %4.9375;
% k_3 = 24.34397227; %parameter(5)*Lce + parameter(6); %17.41*Lce - 2.85;
% k_4 = 18.35536386; %parameter(7)*Lce + parameter(8); %-7.67*Lce + 14.92;
% tau_1 = 0.001895252; %parameter(9); %0.0051;
% tau_2 = 0.015; %parameter(10); % 0.04;
% N = 2.2; %parameter(11)*Lce + parameter(12); %-2.26*Lce + 4.20;
% K = 0.031755997; %parameter(13)*Lce + parameter(14); %-0.044*Lce + 0.080;
% tau_3 = 0.033698122; %parameter(15); %4.475;
% alpha = 1.2;

for i = 237
    MU_No = 69;
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
       
    i
    a = -1;
    b = 1;
    
    r = a + (b-a)*rand(6,1); %0.3
    r = r*0.03;
    r_2 = a + (b-a)*rand(1,1); %0.4
    r_2 = r_2*0.03;
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
    tau_3 = tau_3 + tau_3*r_3;
    
    
    param = [S,C,k_1,k_2,k_3,k_4,tau_1,tau_2,N,K,tau_3,alpha];
    
    simOpt = 0;
    Fs = 10000;
    
    [Data] = MUModel_test_v2(param,Lce,0,MU_type,1,simOpt,Fs);
    MU_No = i;
    cd('/Users/akira/Documents/Github/Twitch-Based-Muscle-Model/Development_Code/Data')
    save(['Data_temp_' num2str(MU_No)],'Data')
    cd(code_folder)
    clear Data
end
