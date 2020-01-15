%==========================================================================
% modelComparison_Fuglevand.m
% Author: Akira Nagamori
% Last update: 6/27/19
% Descriptions:
%   This code is used to generate results_Fuglevand.pdf
%==========================================================================
close all
clear all
clc

%%
condition = 'Model_8_20_CoV_50_Ur_Rec_3';
code_folder = '/Users/akira/Documents/Github/Twitch-Based-Muscle-Model/Code for Figures';
figure_folder = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Figures';

f = 0:0.5:100;

amp_vec = [0.025 0.05 0.1:0.1:1];
data_folder = ['/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Data/withTendon/' condition];
cd(data_folder)
load('mean_Force')
load('std_Force')
load('cov_Force')
load('cov_Force_dt')
load('mean_pxx')
cd(code_folder)
color_code = [37  65 178]/255;

mean_mean_Force = mean(mean_Force);
mean_Force_norm = mean_Force./mean_mean_Force(end)*100;

force_vec = [2 5 15 30 50 70 85 95];
CoV_vec = [6/7*0.5+4.5 2/7*0.5+2.5 5/7*0.5+1 1.5 3.5/7*0.5 + 1.5 6/7*0.5+1 3/7*0.5+1 2/7*0.5+1];


% temp1 = mean(mean_Force_norm);
% temp2 = mean(cov_Force);
% xi = linspace(min(temp1), max(temp1), 50);
% CoV_vec_int = interp1(temp1,temp2,xi,'spline');
% 
% xi_2 = linspace(min(force_vec), max(force_vec), 50);
% CoV_vec_int_2 = interp1(force_vec,CoV_vec,xi_2,'spline');

figure(1)
shadedErrorBar(amp_vec*100,mean(cov_Force_dt),std(cov_Force_dt),'lineprops',{'color',color_code,'LineWidth',2,'markerfacecolor',color_code});
%shadedErrorBar(mean(mean_Force_norm),mean(cov_Force),std(cov_Force),'lineprops',{'color',color_code,'LineWidth',2,'markerfacecolor',color_code});
hold on
plot(force_vec,CoV_vec,'--k','LineWidth',2)
plot(force_vec,CoV_vec,'o','LineWidth',2,'color','k')

% plot(xi,CoV_vec_int,'color',color_code,'LineWidth',2)
% plot(xi_2,CoV_vec_int_2,'--k','LineWidth',2)


figure(1)
%xlabel('Mean Force (%)','FontSize',14)
xlabel('Mean Force (%Maximum Force)','FontSize',14)
ylabel('CoV (%)','FontSize',14)
xlim([0 100])
legend('New model','Moritz et al. 2005','location','northeast')
set(gca,'TickDir','out');
set(gca,'box','off')
% cd (figure_folder)
% saveas(gcf,'activation2CoV_FV_comparison','pdf')
% cd (code_folder)

