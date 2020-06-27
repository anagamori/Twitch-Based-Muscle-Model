close all
clc
clear all

code_folder = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Development_Code';
data_folder = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Model Parameters/Model_11';

MU_type = 'fast';
Lce = 1;

S = 17.7; %parameter(1); %7;
C = 1.8; %parameter(2); %1.025;
k_1 = 65.32; %parameter(3); %14.625;
k_2 = 45.1; %parameter(4); %4.9375;
k_3 = 50.85; %parameter(5)*Lce + parameter(6); %17.41*Lce - 2.85;
k_4 = 34.75; %parameter(7)*Lce + parameter(8); %-7.67*Lce + 14.92;
tau_1 = 0.0023; %parameter(9); %0.0051;
tau_2 = 0.024; %parameter(10); % 0.04;
N = 2.2; %parameter(11)*Lce + parameter(12); %-2.26*Lce + 4.20;
K = 0.039; %parameter(13)*Lce + parameter(14); %-0.044*Lce + 0.080;
tau_3 = 0.07; %parameter(15); %4.475;
alpha = 0.3;

param = [S,C,k_1,k_2,k_3,k_4,tau_1,tau_2,N,K,tau_3,alpha];

simOpt = 1;
Fs = 10000;

[~] = MUModel_test_v2(param,Lce,0,MU_type,1,simOpt,Fs);

%%
%close all
clc
simOpt = 0;
Fs = 2000;

[Data] = MUModel_test_L0(param,Lce,0,MU_type,1,simOpt,Fs);

f_exp = [2 4 6 8 10 12.5 14 16 18 20 22 25 28 30 33.3 40 50 66.6 80 100];
force_exp = [2.5723 7.717 16.399 29.582 44.373 50.482 62.379 63.344 67.203 73.234 77.492 83.923 87.46 88.103 91.318 96.141 97.428 100.32 100 98.392];
[~,loc] = min(abs(force_exp-50));
f_half_exp = f_exp(loc);
fusion_exp = [-0.28169 17.465 36.62 58.31 74.93 85.915 90.141 92.394 94.93 96.056 96.62 97.465 97.465 97.746 98.028 97.746 97.465 97.465 97.746 96.62];

figure(7)
hold on 
plot(f_exp./f_half_exp,fusion_exp/100,'LineWidth',1,'color','k')


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