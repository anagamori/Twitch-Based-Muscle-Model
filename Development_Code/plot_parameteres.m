close all
clc
clear all

code_folder = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Development_Code';
data_folder = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Model Parameters/Model_11';

N_MU = 200;
CT = zeros(N_MU,1);
t2t = zeros(N_MU,1);
FR_half = zeros(N_MU,1);

for i = 1:N_MU

    cd('/Users/akira/Documents/Github/Twitch-Based-Muscle-Model/Development_Code/Data')
    load(['Data_' num2str(i)])
    cd(code_folder)

    CT(i) = Data{2,1};
    t2t(i) = Data{2,5};
    FR_half(i) = Data{2,6};
    
end

figure(1)
plot(CT,t2t,'o')

figure(2)
plot(CT,1./FR_half*1000,'o')

%%
cd('/Users/akira/Documents/Github/Twitch-Based-Muscle-Model/Development_Code/Data')
save('CT','CT')
save('t2t','t2t')
save('FR_half','FR_half')
cd(code_folder)
