close all
clc
clear all

code_folder = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model';
data_folder = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Model Parameters/Model_2/';

N_MU = 300; % number of motor units in a pool
cd(data_folder)
load('index_slow') % index for the largest slow-twitch MU
cd(code_folder)

for trialN = 1:300
    %1:300
    if trialN <= index_slow
        data_folder = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Model Parameters/Model_2/ST';
        MU_type = 'slow';
    else
        data_folder = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Model Parameters/Model_2/FT';
        MU_type = 'fast';
    end
    cd(data_folder)
    load(['MU_' num2str(trialN)])
    cd(code_folder)
    
    [Data] = spikeDrivenMuscleModel_testFunction_fullVersion(parameter,1,0,MU_type,0);
    CT(trialN) = Data{2,1};
    FR_half(trialN) = Data{2,6};
    t2t(trialN) = Data{2,5};
    
end

%%
figure(11)
histogram(t2t)

figure(12)
plot(CT,t2t,'o')
mean(t2t)

CT_vec = CT;
cd('/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Model Parameters/Model_2/')
save('CT_vec','CT_vec')
save('t2t','t2t')
save('FR_half','FR_half')
cd(code_folder)