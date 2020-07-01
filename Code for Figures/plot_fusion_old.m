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

code_folder = '/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Code for Figures';
data_folder = '/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Model Parameters/Model_11';


cd(data_folder)
load('modelParameter')
cd(code_folder)

parameterMatrix = modelParameter.parameterMatrix;

FR_half = zeros(modelParameter.N_MU,52);
Af = zeros(modelParameter.N_MU,52);
fusion = zeros(modelParameter.N_MU,52);
%%
for i = 1:modelParameter.N_MU
    MU_No = i;
    if i <= modelParameter.index_slow
        MU_type = 'slow';
    else
        MU_type = 'fast';
    end
    parameter = parameterMatrix(MU_No,:);

    
    
    [Data] = spikeDrivenMuscleModel_testFunction_fullVersion(parameter,1,0,MU_type,0);
    FR_half(i,:) = Data{2,9};
    Af(i,:) = Data{2,10};
    fusion(i,:) = Data{2,11};
    
    figure(4)
    ax_4 = plot(Data{2,10}*100,Data{2,11}*100,'color',[11,19,43]/255);
    ax_4.Color(4) = 0.5;
    hold on
end

x_scale = 20/2.9;
y_scale =20/2.1;

x = [0.6 2.0 4.6 7 9.5 11.2 12.7 13.5 13.8];
y = [0.25 3 5 7.2 8.8 9.4 9.5 10.1 10.1];

Force = x*x_scale;
Fusion = y*y_scale;

figure(4)
hold on
plot(mean(Af)*100,mean(fusion)*100,'LineWidth',1,'color','k')
plot(Force,Fusion,'LineWidth',2,'color',[252,163,17]/255)
set(gca,'TickDir','out');
set(gca,'box','off')

cd(data_folder)
save('fusion','fusion')
save('Af','Af')
save('FR_half','FR_half')
cd(code_folder)