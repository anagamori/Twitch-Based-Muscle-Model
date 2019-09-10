%==========================================================================
% spikeDrivenMuscleModel_test.m
% Author: Akira Nagamori
% Last update: 5/15/19
% Descriptions
%   A new model driven by spike trains inspired by Williams et al. 1998
%   Generate twitch response and test activation-frequency response with a
%   given parameter set at different muscle lengths
%   Used to generate c-h of Summary_spikeDrivenMuscleModel
%==========================================================================


close all
clc
clear all
%%
code_folder = '/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Code for Figures';
data_folder = '/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Model Parameters/Model_8';

%% 
cd(data_folder)
load('modelParameter')
cd(code_folder)

parameterMatrix = modelParameter.parameterMatrix;

MU_No = 94;
MU_type = 'fast';
parameter = parameterMatrix(MU_No,:);
%%


[Data] = spikeDrivenMuscleModel_testFunction_fullVersion(parameter,1,0,MU_type,1);
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
    [Data] = spikeDrivenMuscleModel_testFunction_fullVersion(parameter,Lce,FR_half,MU_type,1);
    p2p(i+1,:) = Data{2,13};
end

freq = Data{2,9};
fusion = 1-p2p/p2p(5,1);
%%
figure()
plot(freq*FR_half,fusion'*100,'LineWidth',2,'color','b')
xlim([0 100])
xlabel('Frequency (Hz)')
ylabel('Fusion (%)')
%xlim([0 3])
set(gca,'TickDir','out');
set(gca,'box','off')
