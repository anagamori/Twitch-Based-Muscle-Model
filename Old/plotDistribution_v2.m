close all
clc
clear all

code_folder = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model';
data_folder = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Model Parameters/Model_1/FT';

%%
for trialN = 197:300
    %1:300
MU_type = 'fast';

cd(data_folder)
load(['Data_' num2str(trialN)])
cd(code_folder)

cd('/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Model Parameters/Model_1')
load('CT')
cd(code_folder)
%[Data] = spikeDrivenMuscleModel_testFunction_fullVersion(parameter,1,0,MU_type,0);
CT(trialN) = Data{2,1};
FR_half(trialN) = Data{2,6};
t2t(trialN) = Data{2,5};

end

%%
figure(11)
histogram(t2t(197:300))

figure(12)
plot(CT(197:300),t2t(197:300),'o')
mean(t2t(197:300))

