close all
clc
clear all

code_folder = '/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model';
data_folder = '/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Model Parameters/Model_1/ST';

for trialN = 190
    %1:300
MU_type = 'slow';

cd(data_folder)
load(['MU_' num2str(trialN)])
cd(code_folder)

[Data] = spikeDrivenMuscleModel_testFunction_fullVersion(parameter,1,0,MU_type);
FR_half = Data{2,6};
t2t = Data{2,5};
p2p(1,:) = Data{2,13};

for i = 1:4
    if i == 1
        Lce = 0.8;
    elseif i == 2
        Lce = 0.9;
    elseif i == 3
        Lce = 1.1;
    elseif i == 4
        Lce = 1.2;
    end
    [Data] = spikeDrivenMuscleModel_testFunction_fullVersion(parameter,Lce,FR_half,MU_type);
    p2p(i+1,:) = Data{2,13};
end

%close all

end

freq = Data{2,9};
fusion = 1-p2p/p2p(5,1);
plot(freq,fusion','LineWidth',2,'color','b')
xlim([0 3])
