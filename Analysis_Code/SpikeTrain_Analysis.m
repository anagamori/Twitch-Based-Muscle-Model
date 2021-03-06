%==========================================================================
% analysis_muscleModel.m
% Author: Akira Nagamori
% Last update: 3/5/19
% Descriptions:
%   Run muscle model simulation
%==========================================================================
close all
clear all
clc

%%
for k = 1:2
    if k == 1
        condition = 'Model_4_10_CoV_50_Ur_Rec_3';
    else
        condition = 'Model_4_20_CoV_50_Ur_Rec_3';
    end
    data_folder = ['/Volumes/DATA2/New_Model/withTendon/' condition];
    data_folder_git = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Data/withTendon/SpikeData';
    code_folder = '/Users/akira/Documents/Github/Twitch-Based-Muscle-Model/Analysis_Code';
    figure_folder = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Figures';
    
    %%
    for j = 1 %:length(amp_vec)
        if j < 2
            Fs = 10000;
            time = 0:1/Fs:15;
        elseif j >= 2 && j <= 3
            Fs = 15000;
            time = 0:1/Fs:15;
        elseif j >= 4 && j < 7
            Fs = 20000;
            time = 0:1/Fs:15;
        elseif j >= 7
            Fs = 25000;
            time = 0:1/Fs:15;
        end
        
        j
        
        for i = 1
            cd(data_folder)
            load(['Data_' num2str(j) '_' num2str(i)])
            cd(code_folder)
            temp = output.spike_train(1,5*Fs+1:end);
            temp_force = output.force(1,5*Fs+1:end);
            temp_force_2 = output.ForceTendon(5*Fs+1:end);
            
            NumLags = Fs;
            [x,lags] = xcorr(temp,NumLags,'coeff');
            figure(1)
            stem(lags(NumLags+1:end)*1000/Fs,x(NumLags+1:end))
            xlabel('Lags (ms)')
            ylabel('Correlation Coefficient')
            hold on 
            
            [pxx,f] = periodogram(temp-mean(temp),[],0:0.1:100,Fs,'power');
            figure(3)
            plot(f,pxx)
            xlabel('Frequency (Hz)')
            ylabel('Power')
            hold on 
            
            [pxx_force,f] = periodogram(temp_force_2-mean(temp_force_2),[],0:0.1:100,Fs,'power');
            [pxx_force_2,~] = pwelch(temp_force_2-mean(temp_force_2),[],[],0:0.1:100,Fs,'power');
            
            figure(4)
            plot(f,pxx_force./sum(pxx_force)*100)
            hold on 
            plot(f,pxx_force_2./sum(pxx_force_2)*100)
            xlabel('Frequency (Hz)')
            ylabel('Proportion of Total Power (%)')
            hold on 
            
            Data.spike_train = temp;
            Data.unit_force = temp_force;
            Data.force = temp_force_2;
            
            cd(data_folder_git)
            save([condition '_Data'],'Data')
            cd(code_folder)
            
        end
        
    end
end
