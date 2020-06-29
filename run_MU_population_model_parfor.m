%==========================================================================
% run_muscleModel.m
% Author: Akira Nagamori
% Last update: 3/5/19
% Descriptions:
%   Run muscle model simulation
%==========================================================================
close all
clear all
clc


%%
data_folder = '/Volumes/DATA2/PLOS_CB_Data/withTendon/Model_default_v2';
code_folder = '/Users/akira/Documents/Github/Twitch-Based-Muscle-Model';
model_parameter_folder =  '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Development_Code/Data';

%%
cd(model_parameter_folder)
load('modelParameter_v2')
cd(code_folder)

%% Simlulation parameters

amp_vec = [0.025 0.05 0.1:0.1:1];
first = 3;
last = 4;
%parpool(2)
for j = 1:7 %:length(amp_vec)
    j    
    if j <= 7
        Fs = 10000;
        time = 0:1/Fs:15;
    elseif j > 7 && j <= 9
        Fs = 15000;
        time = 0:1/Fs:15;
    elseif j >= 10
        Fs = 20000;
        time = 0:1/Fs:15;
    end
    amp = amp_vec(j);
    input = [zeros(1,1*Fs) amp/2*[0:1/Fs:2] amp*ones(1,length(time)-1*Fs-length(amp*[0:1/Fs:2]))];
    
    Data_temp = cell(1,10);
    tic
    parfor i = first:last
        output_temp = MU_population_model(Fs,time,input,modelParameter,0);
        Data_temp{i} = output_temp;
    end
    toc
    for k = first:last
        output = Data_temp{k};
        cd(data_folder)
        save(['Data_' num2str(j) '_' num2str(k)],'output','-v7.3')
        cd(code_folder)
        figure(1)
        plot(output.ForceTendon)
        hold on
        
        figure(2)
        plot(output.Vce)
        hold on
        Force = output.ForceTendon(end-5*Fs+1:end);
        CoV = std(Force)/mean(Force)
    end
    
    clear output
end

