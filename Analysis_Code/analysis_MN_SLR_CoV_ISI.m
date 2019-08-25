%==========================================================================
% run_muscleModel.m
% Author: Akira Nagamori
% Last update: 8/21/19
% Descriptions:
%   Run muscle model simulation
%==========================================================================
close all
clear all
clc


%%
data_folder = '/Volumes/DATA2/New_Model/SLR/CoV_ISI_test';
code_folder = '/Users/akira/Documents/Github/Twitch-Based-Muscle-Model/Analysis_Code';
model_parameter_folder =  '/Users/akira/Documents/Github/Twitch-Based-Muscle-Model/Model Parameters/Model_7';

%%
cd(model_parameter_folder)
load('modelParameter')
load('parameterMN')
cd(code_folder)

[~,index] = min(modelParameter.U_th);

F0 = modelParameter.F0;
%% 
Fs = 30000;
%%
amp_vec = [0.01:0.01:1];
lastTrial = 30;

%%
spike_mat = zeros(lastTrial,10*Fs+1);

for j = 1:lastTrial %length(amp_vec)
    j
    Fs = 30000;
    time = 0:1/Fs:15;
    i = 1;
    
    tic
    cd(data_folder)
    load(['Data_' num2str(j) '_' num2str(i)])
    cd(code_folder)
    toc
    
    temp = output.ForceTendon(5*Fs+1:end);
    [pxx,f] = pwelch(temp-mean(temp),[],[],0:0.1:100,Fs,'power');
    figure(11)
    plot(f,pxx,'LineWidth',2)
    hold on
    xlim([0 30])
    xlabel('Frequency (Hz)','FontSize',14)
    ylabel('Power (N^2)','FontSize',14)
    
    mean_Force(j) = mean(temp);
    Force_CoV(j) = std(temp)/mean(temp)*100;
    
    temp_spikeTrain = output.spike_train(index,5*Fs+1:end);
    spike_mat(j,:) = temp_spikeTrain;
    
    clear output
end

figure(1)
plot(amp_vec(1:lastTrial)*100,mean_Force./F0*100)
xlabel('Activation (%Maximum)','FontSize',14)
ylabel('%MVC (%)','FontSize',14)

figure(2)
plot(amp_vec(1:lastTrial)*100,Force_CoV)
xlabel('Activation (%Maximum)','FontSize',14)
ylabel('CoV for Force (%)','FontSize',14)

%%
for i = 1:lastTrial
    spike_time = find(spike_mat(i,:));
    ISI = diff(spike_time)/(Fs/1000);
    FR_vec(i) = nanmean(1./ISI*1000);
    FR_sd(i) = nanstd(1./ISI*1000);
    FR_cov(i) = FR_sd(i)/FR_vec(i)*100;
end
figure(4)
plot(amp_vec(1:lastTrial)*100,FR_cov)
xlabel('Activation (%Maximum)','FontSize',14)
ylabel('CoV for Force (%)','FontSize',14)

