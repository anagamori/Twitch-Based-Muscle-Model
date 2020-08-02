%==========================================================================
% run_muscleModel.m
% Author: Akira Nagamori
% Last update: 3/5/19
% Descriptions:
%   Run muscle model simulation
%==========================================================================
%close all
clear all
clc


%%
data_folder = 'D:\MU_Population_Model\Old_model_v3';
code_folder = 'C:\Users\anaga\OneDrive\Documents\GitHub\Twitch-Based-Muscle-Model';
model_parameter_folder =  'C:\Users\anaga\OneDrive\Documents\GitHub\Twitch-Based-Muscle-Model\Development_Code\Data';

%%
cd(model_parameter_folder)
load('modelParameter_v2_length_v2')
cd(code_folder)

%% Simlulation parameters

amp_vec = [0.025 0.05 0.1:0.1:1];
%amp_vec = [0.106 0.29 0.62 0.78 0.93];
trial_vec = [7 10];
for j = 12 %1:length(amp_vec)
    j
    %      Fs = 20000;
    %         time = 0:1/Fs:15;
    if j <= 7
        Fs = 10000;
        time = 0:1/Fs:5;
    elseif j > 7 && j <= 9
        Fs = 15000;
        time = 0:1/Fs:5;
    elseif j >= 10
        Fs = 30000;
        time = 0:1/Fs:10;
    end
    
    amp = amp_vec(j);
    input = [zeros(1,1*Fs) amp/2*[0:1/Fs:2] amp*ones(1,length(time)-1*Fs-length(amp*[0:1/Fs:2]))];
    %%
    for i = 1
        i
        tic
        output = MU_population_model(Fs,time,input,modelParameter,1);
        toc
        cd(data_folder)
        save(['Data_' num2str(j) '_' num2str(i)],'output','-v7.3')
        cd(code_folder)
        Force = output.ForceTendon(end-1*Fs+1:end);
        mean(Force(end-1*Fs+1:end))/444.5844
        CoV = std(Force)/mean(Force)
        %        clear output
        mean_FR = zeros(1,200);
        CoV_FR = zeros(1,200);
        
        for n = 1:200
            spike_time = find(output.spike_train(n,end-1*Fs+1:end));
            ISI = diff(spike_time)/(Fs/1000);
            mean_FR(n) = mean(1./ISI*1000);
            CoV_FR(n) = std(ISI)/mean(ISI)*100;
        end
    end
    
end

figure(11)
plot(mean_FR,'o')
hold on

figure(12)
plot(CoV_FR,'o')
hold on 