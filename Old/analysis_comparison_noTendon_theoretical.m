%==========================================================================
% analysis_comparison.m
% Author: Akira Nagamori
% Last update: 3/11/19
% Descriptions:
%==========================================================================
close all
clear all
clc

%%
condition = '10_CoV_50_Ur_Rec_2';
code_folder = '/Users/akiranagamori/Documents/Github/Twitch-Based-Muscle-Model';
figure_folder = '/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Figures';


amp_vec = 0.1:0.1:1;
mean_Force = zeros(10,length(amp_vec));
std_Force = zeros(10,length(amp_vec));
cov_Force = zeros(10,length(amp_vec));
pxx = zeros(10,201);
mean_pxx = zeros(length(amp_vec),201);
f = 0:0.5:100;

for i = 1
    if i == 1
        Fs = 1000;
        time =0:1/Fs:15;
        data_folder = ['/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Data/noTendon/' condition];
        cd(data_folder)
        load('mean_Force')
        load('std_Force')
        load('cov_Force')
        load('mean_pxx')
        cd(code_folder)
    elseif i == 2
        Fs = 10000;
        time =0:1/Fs:15;
        data_folder = ['/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Data/withTendon/' condition];
        cd(data_folder)
        load('mean_Force')
        load('std_Force')
        load('cov_Force')
        load('mean_pxx')
        cd(code_folder)
    end
    mean_mean_Force = mean(mean_Force);
    figure(1)
    plot([0 amp_vec],[0 mean(mean_Force)]./mean_mean_Force(end),'LineWidth',2)
    %plot([0 amp_vec],[0 mean(mean_Force)],'LineWidth',2)
    xlabel('Activation (%)','FontSize',14)
    ylabel('Force (%MVC)','FontSize',14)
    set(gca,'TickDir','out');
    set(gca,'box','off')
    hold on
    
    figure(2)
    plot([0 mean(mean_Force)./mean_mean_Force(end)*100],[0 mean(std_Force/mean_mean_Force(end)*100)],'LineWidth',2);
    hold on
    plot(0:1:100,[0:1:100]*(3/100))
    
    %errorbar(amp_vec,mean(std_Force),std(std_Force),'LineWidth',2);
    xlabel('Mean Force (%)','FontSize',14)
    ylabel('SD (N)','FontSize',14)
    hold on
    yticks([0.05 0.1 0.15 0.2 0.25])
    set(gca,'TickDir','out');
    set(gca,'box','off')

    figure(3)
    errorbar(mean(mean_Force)./mean_mean_Force(end),mean(cov_Force),std(cov_Force),'LineWidth',2);
    %errorbar(amp_vec,mean(cov_Force),std(cov_Force),'LineWidth',2);
    xlabel('Mean Force (%)','FontSize',14)
    ylabel('CoV (%)','FontSize',14)
    hold on
    set(gca,'TickDir','out');
    set(gca,'box','off')
    
    
end

x = log([0 mean(mean_Force)./mean_mean_Force(end)*100])';
y = log([0 mean(std_Force/mean_mean_Force(end)*100)])';