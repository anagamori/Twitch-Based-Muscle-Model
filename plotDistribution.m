close all
clc
clear all

code_folder = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model';
data_folder = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Model Parameters/Model_1/ST';

for trialN = 1:100
    %1:300
MU_type = 'slow';

cd(data_folder)
load(['MU_' num2str(trialN)])
cd(code_folder)

[Data] = spikeDrivenMuscleModel_testFunction_fullVersion(parameter,1,0,MU_type);
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

