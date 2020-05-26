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
data_folder = '/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Model Parameters/Model_11';

%% 
cd(data_folder)
load('modelParameter')
cd(code_folder)

parameterMatrix = modelParameter.parameterMatrix;

MU_No = 30;
MU_type = 'slow';
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
%% Data from McNulty et al. (2000)
f_exp = [2 4 6 8 10 12.5 14 16 18 20 22 25 28 30 33.3 40 50 66.6 80 100];
force_exp = [2.5723 7.717 16.399 29.582 44.373 50.482 62.379 63.344 67.203 73.234 77.492 83.923 87.46 88.103 91.318 96.141 97.428 100.32 100 98.392];
[~,loc] = min(abs(force_exp-50));
f_half_exp = f_exp(loc);
fusion_exp = [-0.28169 17.465 36.62 58.31 74.93 85.915 90.141 92.394 94.93 96.056 96.62 97.465 97.465 97.746 98.028 97.746 97.465 97.465 97.746 96.62];

%% Data from Macefield et al. (1993)
% f_half_exp =  10.741;
% f_exp = [2 5 8 10 15 20 30 50 80 100];
% force_exp = [4.1775 14.36 31.854 44.909 64.752 76.762 86.945 91.906 92.689 94.256];
% fusion_exp = [3.7123 27.842 48.956 67.981 83.527 88.167 88.631 95.128 95.592 95.592];
%%
figure()
plot(freq*FR_half,fusion'*100,'LineWidth',2,'color','b')
xlim([0 100])
xlabel('Frequency (Hz)')
ylabel('Fusion (%)')
%xlim([0 3])
set(gca,'TickDir','out');
set(gca,'box','off')

%%
figure()
plot(freq,(1-p2p(4,:)./p2p(4,1))*100,'LineWidth',1,'color','b')
hold on 
plot(f_exp./f_half_exp,fusion_exp,'LineWidth',1,'color','k')
xlim([0 3])
xlabel('Frequency (Hz)')
ylabel('Fusion (%)')
%xlim([0 3])
set(gca,'TickDir','out');
set(gca,'box','off')