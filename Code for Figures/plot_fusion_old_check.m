
close all
clc
clear all

code_folder = '/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Code for Figures';
data_folder = '/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Model Parameters/Model_11';


cd(data_folder)
load('modelParameter')
load('fusion')
load('Af')
load('FR_half')
cd(code_folder)

parameterMatrix = modelParameter.parameterMatrix;

%%
idx_nan_temp = isnan(Af);
idx_nan = find(idx_nan_temp);
[row,col] = ind2sub(size(Af),idx_nan);

%%
for i = 136
    MU_No = i;
    if i <= modelParameter.index_slow
        MU_type = 'slow';
    else
        MU_type = 'fast';
    end
    parameter = parameterMatrix(MU_No,:);


    [Data] = spikeDrivenMuscleModel_testFunction_fullVersion(parameter,1,0,MU_type,1);
    
    FR_half(i,:) = Data{2,9};
    Af(i,:) = Data{2,10};
    fusion(i,:) = Data{2,11};
    
    figure(4)
    ax_4 = plot(Data{2,10}*100,Data{2,11}*100,'color',[11,19,43]/255);
    ax_4.Color(4) = 0.5;
    hold on
end