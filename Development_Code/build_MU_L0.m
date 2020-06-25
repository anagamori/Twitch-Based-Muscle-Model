close all
clc
clear all

code_folder = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Development_Code';
data_folder = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Model Parameters/Model_11';

MU_type = 'fast';
Lce = 1;

S = 19; %parameter(1); %7;
C = 1.8; %parameter(2); %1.025;
k_1 = 75.8; %parameter(3); %14.625;
k_2 = 72.6; %parameter(4); %4.9375;
k_3 = 65.2; %parameter(5)*Lce + parameter(6); %17.41*Lce - 2.85;
k_4 = 63; %parameter(7)*Lce + parameter(8); %-7.67*Lce + 14.92;
tau_1 = 0.002; %parameter(9); %0.0051;
tau_2 = 0.023; %parameter(10); % 0.04;
N = 2.2; %parameter(11)*Lce + parameter(12); %-2.26*Lce + 4.20;
K = 0.049; %parameter(13)*Lce + parameter(14); %-0.044*Lce + 0.080;
tau_3 = 0.035; %parameter(15); %4.475;
alpha = 1.2;

param = [S,C,k_1,k_2,k_3,k_4,tau_1,tau_2,N,K,tau_3,alpha];

simOpt = 1;
Fs = 10000;

[~] = MUModel_test_v2(param,Lce,0,MU_type,1,simOpt,Fs);

%%
close all
clc
simOpt = 0;
Fs = 2000;

[Data] = MUModel_test_v2(param,Lce,0,MU_type,1,simOpt,Fs);

%%
close all
clc
simOpt = 0;
Fs = 10000;

[Data] = MUModel_test_v2(param,Lce,0,MU_type,1,simOpt,Fs);

%%
MU_No = 56;
cd('/Users/akira/Documents/Github/Twitch-Based-Muscle-Model/Development_Code/Data')
save(['Data_' num2str(MU_No)],'Data')
cd(code_folder)
% cd('/Users/akira/Documents/Github/Twitch-Based-Muscle-Model/Development_Code/Data')
% save('seed_ST_1','parameter')
% cd(code_folder)