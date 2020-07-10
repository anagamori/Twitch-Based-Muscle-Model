%==========================================================================
% comparison_CoV_ISI.m
% Author: Akira Nagamori
% Last update: 7/18/19
% Descriptions:
%   This code is used to generate fig XX in the paper
%==========================================================================
close all
clear all
clc

%%
code_folder = '/Users/akira/Documents/Github/Twitch-Based-Muscle-Model/Code for Figures';
figure_folder = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Figures';

Fs = 10000;
T = 0.07;
t_twitch = 0:1/Fs:1;
twitch = t_twitch./T.*exp(1-t_twitch./T);

T_2 = 0.06;
t_twitch = 0:1/Fs:1;
twitch_2 = t_twitch./T_2.*exp(1-t_twitch./T_2);

[pxx_twitch,f] = pwelch(twitch,[],[],0:0.1:100,Fs,'power');
pxx_twitch_norm = pxx_twitch./sum(pxx_twitch);

unit = 4;
CoV = zeros(3,10);
tic
for i = 1:3
    pxx_force = zeros(10,length(0:0.1:100));
    CoV_force = zeros(10,1);
    if i == 1
        data_folder = '/Volumes/DATA2/PLOS_CB_Data/withTendon/Individual_Unit/CoV_0';
        color_code = [0 0 0];
    elseif i == 2
        data_folder = '/Volumes/DATA2/PLOS_CB_Data/withTendon/Individual_Unit/CoV_10';
        color_code = [230 57 70]/255;
    elseif i == 3
        data_folder = '/Volumes/DATA2/PLOS_CB_Data/withTendon/Individual_Unit/CoV_20';
        color_code = [37  65 178]/255;
    end
    for j = 1:10
        cd(data_folder)
        load(['Data_' num2str(2) '_' num2str(j)])
        cd(code_folder)
        Force = output.ForceTendon(5*Fs+1:end);
        [pxx_force(j,:),~] = pwelch(Force-mean(Force),gausswin(5*Fs),0.9*5*Fs,0:0.1:100,Fs);
        CoV_force(j) = std(Force)/mean(Force)*100;
    end
    CoV(i,:) = CoV_force';
    
    %%
    figure(1)
    plot(f,mean(pxx_force),'LineWidth',1,'color',color_code)
    hold on
    xlim([0 30])
    xlabel('Frequency (Hz)','FontSize',10)
    ylabel('Power (N^2/Hz)','FontSize',10)
    set(gca,'TickDir','out');
    set(gca,'box','off')
end
toc
figure(2)
boxplot(CoV')
xlabel('Frequency (Hz)','FontSize',10)
ylabel('Power (N^2/Hz)','FontSize',10)
set(gca,'TickDir','out');
set(gca,'box','off')
