%==========================================================================
% analysis_comparison.m
% Author: Akira Nagamori
% Last update: 4/11/19
% Descriptions:
%   This code is used to generate fig XX in the paper
%==========================================================================
close all
clear all
clc

%%
condition = 'Model_4_10_CoV_50_Ur_Rec_3';
code_folder = '/Users/akiranagamori/Documents/Github/Twitch-Based-Muscle-Model';
figure_folder = '/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Figures';

amp_vec = [0.05 0.1:0.1:1];
mean_Force = zeros(10,length(amp_vec));
std_Force = zeros(10,length(amp_vec));
cov_Force = zeros(10,length(amp_vec));
pxx = zeros(10,201);
mean_pxx = zeros(length(amp_vec),201);
f = 0:0.5:100;

CoV = [];
for i = 1:2
    if i == 1
        Fs = 1000;
        time =0:1/Fs:15;
        data_folder = ['/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Data/SLR/Ia_gain_10000' ];
        cd(data_folder)
        load('mean_Force')
        load('std_Force')
        load('cov_Force')
        load('mean_pxx')
        cd(code_folder)
        color_code = [217 4 41]/255;
        cov_2 = cov_Force;
    elseif i == 2
        Fs = 10000;
        time =0:1/Fs:15;
        data_folder = ['/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Data/SLR/Ia_gain_2000'];
        cd(data_folder)
        load('mean_Force')
        load('std_Force')
        load('cov_Force')
        load('mean_pxx')
        cd(code_folder)
        color_code = [37  65 178]/255;
        cov_2 = cov_Force;
    end
    mean_mean_Force = mean(mean_Force);
    CoV = [CoV cov_Force(:,2)];
    
    %%
    figure(4)
    plot(f,mean_pxx(2,:),'LineWidth',2,'color',color_code)
    hold on
    xlim([0 50])
    xlabel('Frequency (Hz)','FontSize',10)
    ylabel('Power (N^2)','FontSize',10)
    set(gca,'TickDir','out');
    set(gca,'box','off')
    ax = gca;
    ax.FontSize = 6;
    
end

%%
figure(1)
boxplot(CoV)
xlabel('Ia Gain','FontSize',14)
ylabel('CoV (%)','FontSize',14)
set(gca,'TickDir','out');
set(gca,'box','off')
legend('Without Tendon','With Tendon','With Tendon & no FV','Shorter Tendon','location','northwest')

