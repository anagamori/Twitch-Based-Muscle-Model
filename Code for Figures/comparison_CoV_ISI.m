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
CoV = zeros(10,2);
for i = 1:2
    if i == 1
        condition = 'Model_4_10_CoV_50_Ur_Rec_3';
        data_folder = ['/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Data/withTendon/SpikeData'];
        cd(data_folder)
        load([condition '_Data'])
        cd(code_folder)
        
        data_folder = ['/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Data/withTendon/' condition];
        cd(data_folder)
        load('cov_Force')
        load('mean_pxx')
        cd(code_folder)
        
        color_code = [230 57 70]/255;
        
    elseif i == 2
        condition = 'Model_4_20_CoV_50_Ur_Rec_3';
        data_folder = ['/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Data/withTendon/SpikeData'];
        cd(data_folder)
        load([condition '_Data'])
        cd(code_folder)
        
        data_folder = ['/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Data/withTendon/' condition];
        cd(data_folder)
        load('cov_Force')
        load('mean_pxx')
        cd(code_folder)
        
        color_code = [37  65 178]/255;
    end
    
    %%
    CoV(:,i) = cov_Force(:,2);
    %%
    [pxx_spike,f] = pwelch(Data.spike_train-mean(Data.spike_train),[],[],0:0.5:100,Fs,'power');
    [pxx_unit_force,~] = pwelch(Data.unit_force-mean(Data.unit_force),[],[],0:0.5:100,Fs,'power');
    [pxx_force,~] = pwelch(Data.force-mean(Data.force),[],[],0:0.5:100,Fs,'power');
    %%
    figure(4)
    subplot(2,2,1);
    plot(f,pxx_spike./sum(pxx_spike)*100,'LineWidth',1,'color',color_code)
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
    
    subplot(2,2,2);
    plot(f,pxx_unit_force./sum(pxx_unit_force)*100,'LineWidth',1,'color',color_code)
    hold on
    xlim([0 50])
    yticks(0:10:50)
    xlabel('Frequency (Hz)','FontSize',10)
    ylabel('Proportion of Total Power (%)','FontSize',10)
    set(gca,'TickDir','out');
    set(gca,'box','off')
    ax = gca;
    ax.FontSize = 6;
    hold on
    
    subplot(2,2,3);
    plot(f,mean_pxx(2,:)./sum(mean_pxx(2,:))*100,'LineWidth',1,'color',color_code)
    hold on
    xlim([0 50])
    xlabel('Frequency (Hz)','FontSize',10)
    ylabel('Proportion of Total Power (%)','FontSize',10)
    set(gca,'TickDir','out');
    set(gca,'box','off')
    ax = gca;
    ax.FontSize = 6;
    
%     subplot(2,2,4);
%     plot(f,mean_pxx(8,:),'LineWidth',2,'color',color_code)
%     xlim([0 50])
%     hold on
%     xlabel('Frequency (Hz)','FontSize',10)
%     ylabel('Power (N^2)','FontSize',10)
%     set(gca,'TickDir','out');
%     set(gca,'box','off')
%     ax = gca;
%     ax.FontSize = 6;
    
    
end

figure(4)
subplot(2,2,4)
boxplot(CoV,'Labels',{'CoV ISI = 10%','CoV ISI = 20%'})
ylabel('CoV for Force (%)','FontSize',10)
yticks(0.6:0.2:1.4)
set(gca,'TickDir','out');
set(gca,'box','off')
ax = gca;
ax.FontSize = 6;

figure(4)
fig = gcf;
%linkaxes([ax1,ax2,ax3,ax4],'y')
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6.92 6.92];
