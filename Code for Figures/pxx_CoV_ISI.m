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

unit = 10;
CoV = zeros(10,2);
for i = 1:3
    if i == 1
  
        data_folder = '/Volumes/DATA2/New_Model/withTendon/Model_8_no_CoV_50_Ur_Rec_3';
        cd(data_folder)
        load(['Data_' num2str(0) '_' num2str(1)])
        cd(code_folder)      
        color_code = [0 0 0];       
    elseif i == 2
     
        data_folder = '/Volumes/DATA2/New_Model/withTendon/Model_8_10_CoV_50_Ur_Rec_3';
        cd(data_folder)
        load(['Data_' num2str(0) '_' num2str(1)])
        cd(code_folder)
        color_code = [230 57 70]/255;              
    elseif i == 3
     
        data_folder = '/Volumes/DATA2/New_Model/withTendon/Model_8_20_CoV_50_Ur_Rec_3';
        cd(data_folder)
        load(['Data_' num2str(0) '_' num2str(1)])
        cd(code_folder)              
        color_code = [37  65 178]/255;
    end
    spike_train = output.spike_train(unit,:);
    Force = output.ForceTendon;
    force = output.force(unit,:);
    
    force_temp = conv(spike_train,twitch);
    force_temp = force_temp(5*Fs+1:length(spike_train));
    cov(i) = std(force_temp)/mean(force_temp)*100;
    
    force_temp_2 = conv(spike_train,twitch_2);
    force_temp_2 = force_temp_2(5*Fs+1:length(spike_train));
    cov_2(i) = std(force_temp_2)/mean(force_temp_2)*100;
    
    figure(11)
    plot(force_temp)
    hold on
    plot(force_temp_2)
    %%
    [pxx_spike,~] = pwelch(spike_train(5*Fs+1:end)-mean(spike_train(5*Fs+1:end)),[],[],0:0.1:100,Fs);
    pxx_spike_norm = pxx_spike./sum(pxx_spike);
    [pxx_unit_force,~] = pwelch(force(5*Fs+1:end)-mean(force(5*Fs+1:end)),[],[],0:0.1:100,Fs);
    [pxx_force,~] = pwelch(Force(5*Fs+1:end)-mean(Force(5*Fs+1:end)),[],[],0:0.1:100,Fs);
    
    [pxx_unit_force_conv,~] = pwelch(force_temp-mean(force_temp),[],[],0:0.1:100,Fs);
    [pxx_unit_force_conv_2,~] = pwelch(force_temp_2-mean(force_temp_2),[],[],0:0.1:100,Fs);
    %%
    figure(1)
    plot(f,pxx_spike,'LineWidth',1,'color',color_code)
    hold on
    xlim([0 30])
    yticks(0:1:4)
    title('Unit Force')
    xlabel('Frequency (Hz)','FontSize',10)
    ylabel('Proportion of Total Power (%)','FontSize',10)
    set(gca,'TickDir','out');
    set(gca,'box','off')
    
    figure(2)
    plot(f,pxx_spike,'LineWidth',1,'color',color_code)
    hold on
    xlim([0 10])
    %yticks(0:1:4)
    title('Unit Force')
    xlabel('Frequency (Hz)','FontSize',10)
    ylabel('Proportion of Total Power (%)','FontSize',10)
    set(gca,'TickDir','out');
    set(gca,'box','off')
%     ax = gca;
%     ax.FontSize = 6;      
    
    pxx_conv = pxx_spike_norm.*pxx_twitch_norm;
    figure(3)
    plot(f,pxx_conv./sum(pxx_conv)*100,'LineWidth',1,'color',color_code)
    hold on
    xlim([0 30])
    title('Convolution with a fixed twitch profile')
    xlabel('Frequency (Hz)','FontSize',10)
    ylabel('Proportion of Total Power (%)','FontSize',10)
    set(gca,'TickDir','out');
    set(gca,'box','off')
    
    %./sum(pxx_unit_force)*100
    figure(4)
    plot(f,pxx_unit_force,'LineWidth',1,'color',color_code)
    hold on
    xlim([0 30])
    %yticks(0:10:50)
    title('Spke Train')
    xlabel('Frequency (Hz)','FontSize',10)
    ylabel('Proportion of Total Power (%)','FontSize',10)
    set(gca,'TickDir','out');
    set(gca,'box','off')
%     ax = gca;
%     ax.FontSize = 6;
    hold on
    
    figure(5)
    plot(f,pxx_force,'LineWidth',1,'color',color_code)
    hold on
    xlim([0 30])
    %yticks(0:10:50)
    title('Spke Train')
    xlabel('Frequency (Hz)','FontSize',10)
    ylabel('Proportion of Total Power (%)','FontSize',10)
    set(gca,'TickDir','out');
    set(gca,'box','off')
end
figure(1)
legend('CoV ISI = 0%','CoV ISI = 10%','CoV ISI = 20%')

figure(2)
legend('CoV ISI = 0%','CoV ISI = 10%','CoV ISI = 20%')

