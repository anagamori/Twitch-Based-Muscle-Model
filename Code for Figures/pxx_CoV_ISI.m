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
code_folder = '/Users/akiranagamori/Documents/Github/Twitch-Based-Muscle-Model/Code for Figures';
figure_folder = '/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Figures';

Fs = 10000;
T = 0.09;
t_twitch = 0:1/Fs:1;
twitch = t_twitch./T.*exp(1-t_twitch./T);
[pxx_twitch,f] = pwelch(twitch,[],[],0:0.1:100,Fs,'power');
pxx_twitch_norm = pxx_twitch./sum(pxx_twitch);

CoV = zeros(10,2);
for i = 1:2
    if i == 1
        condition = 'Model_4_10_CoV_50_Ur_Rec_3';
        data_folder = ['/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Data/withTendon/SpikeData'];
        cd(data_folder)
        load([condition '_Data'])
        cd(code_folder)
        
        color_code = [230 57 70]/255;
        
    elseif i == 2
        condition = 'Model_4_20_CoV_50_Ur_Rec_3';
        data_folder = ['/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Data/withTendon/SpikeData'];
        cd(data_folder)
        load([condition '_Data'])
        cd(code_folder)
              
        color_code = [37  65 178]/255;
    end
    

    %%
    [pxx_spike,~] = pwelch(Data.spike_train-mean(Data.spike_train),[],[],0:0.1:100,Fs,'power');
    pxx_spike_norm = pxx_spike./sum(pxx_spike);
    [pxx_unit_force,~] = pwelch(Data.unit_force-mean(Data.unit_force),[],[],0:0.1:100,Fs,'power');
    [pxx_force,~] = pwelch(Data.force-mean(Data.force),[],[],0:0.5:100,Fs,'power');
    %%
    figure(1)
    plot(f,pxx_spike_norm*100,'LineWidth',1,'color',color_code)
    hold on
    xlim([0 50])
    yticks(0:1:4)
    xlabel('Frequency (Hz)','FontSize',10)
    ylabel('Proportion of Total Power (%)','FontSize',10)
    set(gca,'TickDir','out');
    set(gca,'box','off')
    ax = gca;
    ax.FontSize = 6;
    legend('CoV ISI = 10%','CoV ISI = 20%')
    
    figure(2)
    plot(f,pxx_twitch_norm*100,'LineWidth',1,'color',color_code)
    hold on
    xlim([0 50])
    yticks(0:1:4)
    xlabel('Frequency (Hz)','FontSize',10)
    ylabel('Proportion of Total Power (%)','FontSize',10)
    set(gca,'TickDir','out');
    set(gca,'box','off')
    ax = gca;
    ax.FontSize = 6;
    legend('CoV ISI = 10%','CoV ISI = 20%')
    
    pxx_conv = pxx_spike_norm.*pxx_twitch_norm;
    figure(3)
    plot(f,pxx_conv./sum(pxx_conv)*100,'LineWidth',1,'color',color_code)
    hold on
    xlim([0 50])
    %yticks(0:10:50)
    xlabel('Frequency (Hz)','FontSize',10)
    ylabel('Proportion of Total Power (%)','FontSize',10)
    set(gca,'TickDir','out');
    set(gca,'box','off')
    ax = gca;
    ax.FontSize = 6;
    hold on
    
    figure(4)
    plot(f,pxx_unit_force./sum(pxx_unit_force)*100,'LineWidth',1,'color',color_code)
    hold on
    xlim([0 50])
    %yticks(0:10:50)
    xlabel('Frequency (Hz)','FontSize',10)
    ylabel('Proportion of Total Power (%)','FontSize',10)
    set(gca,'TickDir','out');
    set(gca,'box','off')
    ax = gca;
    ax.FontSize = 6;
    hold on
    
    
end

