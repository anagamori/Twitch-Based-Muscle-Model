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
data_folder = '/Volumes/DATA2/PLOS_CB_Data/withTendon/Model_N_400';
code_folder = '/Users/akira/Documents/Github/Twitch-Based-Muscle-Model';
model_parameter_folder =  '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Development_Code/Data';

%%
cd(model_parameter_folder)
load('modelParameter_N_400')
cd(code_folder)

%% Simlulation parameters

amp_vec = [0.025 0.05 0.1:0.1:1];
trial_vec = [7 10];
for j = 8:length(amp_vec)
    j
    %     if j <= 1+2
    %         Fs = 10000;
    %         time = 0:1/Fs:5;
    %     elseif j >= 2+2 && j <= 4+2
    %         Fs = 10000;
    %         time = 0:1/Fs:5;
    %     elseif j > 4+2 && j < 7+2
    %         Fs = 20000;
    %         time = 0:1/Fs:5;
    %     elseif j >= 7+2 && j <= 8+2
    %         Fs = 15000;
    %         time = 0:1/Fs:5;
    %     elseif j >= 9+2
    %         Fs = 20000;
    %         time = 0:1/Fs:5;
    %     end
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
    %     Fs = 10000;
    %         time = 0:1/Fs:15;
    amp = amp_vec(j);
    input = [zeros(1,1*Fs) amp/2*[0:1/Fs:2] amp*ones(1,length(time)-1*Fs-length(amp*[0:1/Fs:2]))];
    %input = [zeros(1,1*Fs) amp/2*[0:1/Fs:2] amp*ones(1,2*Fs) amp-amp/50*[1:10*Fs]/Fs];
    %%
    for i = 2:3
        i
        tic
        output = MU_population_model(Fs,time,input,modelParameter,1);
        toc
        cd(data_folder)
        save(['Data_' num2str(j) '_' num2str(i)],'output','-v7.3')
        cd(code_folder)
        Force = output.ForceTendon(end-5*Fs+1:end);
        CoV = std(Force)/mean(Force)
        clear output
    end
    %     if j == 1
    %         for i = 8:10
    %             i
    %             tic
    %             output = MU_population_model(Fs,time,input,modelParameter,1);
    %             toc
    %             cd(data_folder)
    %             save(['Data_' num2str(j) '_' num2str(i)],'output','-v7.3')
    %             cd(code_folder)
    %         end
    %     elseif j == 10
    %         for i = 8:10
    %             i
    %             tic
    %             output = MU_population_model(Fs,time,input,modelParameter,1);
    %             toc
    %             cd(data_folder)
    %             save(['Data_' num2str(j) '_' num2str(i)],'output','-v7.3')
    %             cd(code_folder)
    %
    %         end
    %     else
    %         for i = 3
    %             i
    %             tic
    %             output = MU_population_model(Fs,time,input,modelParameter,1);
    %             toc
    %             cd(data_folder)
    %             save(['Data_' num2str(j) '_' num2str(i)],'output','-v7.3')
    %             cd(code_folder)
    %
    %
    %         end
    %     end
    
    %     figure(2)
    %     plot(output.Vce)
    %     hold on
    
    
end

