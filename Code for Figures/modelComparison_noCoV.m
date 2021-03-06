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
condition = 'Model_8_no_CoV_50_Ur_Rec_3';
code_folder = '/Users/akiranagamori/Documents/Github/Twitch-Based-Muscle-Model/Code for Figures';
figure_folder = '/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Figures';

f = 0:0.5:100;

amp_vec = [0.05 0.1:0.1:1];
data_folder = ['/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Data/withTendon/' condition];
cd(data_folder)
load('mean_Force')
load('std_Force')
load('cov_Force')
load('mean_pxx')
cd(code_folder)
color_code = [37  65 178]/255;

mean_mean_Force = mean(mean_Force);
mean_Force_norm = mean_Force./mean_mean_Force(end)*100;

force_vec = [2 5 15 30 50 70 85 95];
CoV_vec = [6/7*0.5+4.5 2/7*0.5+2.5 5/7*0.5+1 1.5 3.5/7*0.5 + 1.5 6/7*0.5+1 3/7*0.5+1 2/7*0.5+1];


figure(1)
%shadedErrorBar(amp_vec*100,mean(cov_Force),std(cov_Force),'lineprops',{'color',color_code,'LineWidth',2,'markerfacecolor',color_code});
%shadedErrorBar(mean(mean_Force_norm),mean(cov_Force),std(cov_Force),'lineprops',{'color',color_code,'LineWidth',2,'markerfacecolor',color_code});
plot(amp_vec*100,cov_Force,'color',color_code,'LineWidth',2)
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

