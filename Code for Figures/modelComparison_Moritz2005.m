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
condition = 'Model_8_var_CoV_50_Ur_Rec_3';
code_folder = '/Users/akira/Documents/Github/Twitch-Based-Muscle-Model/Code for Figures';
figure_folder = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Figures';

f = 0:0.5:100;

amp_vec = [0.025 0.05 0.1:0.1:1];
data_folder = ['/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Data/withTendon/' condition];
cd(data_folder)
load('mean_Force')
load('std_Force')
load('std_Force_dt')
load('cov_Force')
load('cov_Force_dt')
load('mean_pxx')
cd(code_folder)
color_code = [37  65 178]/255;

mean_mean_Force = mean(mean_Force);
mean_Force_norm = mean_Force./mean_mean_Force(end)*100;
max_Force = mean_mean_Force(end);

force_vec = [2 5 15 30 50 70 85 95];
std_vec = [6/14*0.25 8/14*0.25 12/14*0.25 11/14*0.25+0.25 10/14*0.25+0.75 2/14*0.25+1 4.5/14*0.25+1 10/14*0.25+1];
CoV_vec = [6/7*0.5+4.5 2/7*0.5+2.5 5/7*0.5+1 1.5 3.5/7*0.5 + 1.5 6/7*0.5+1 3/7*0.5+1 2/7*0.5+1];

std_Force_dt = std_Force_dt/max_Force*100;
 
a = std_vec(1)/force_vec(1);
a_todorov = 0.1289;
x = [0 amp_vec]*100;


figure(1)
shadedErrorBar(mean(mean_Force_norm),mean(cov_Force_dt),std(cov_Force_dt),'lineprops',{'color',color_code,'LineWidth',2,'markerfacecolor',color_code});
hold on
plot(force_vec,CoV_vec,'--k','LineWidth',2)
plot(force_vec,CoV_vec,'o','LineWidth',2,'color','k')
plot([x(1) x(end)],[0.02*100 0.02*100],'LineWidth',2,'Color',[77 172 38]/255)
%plot([x(1) x(end)],[a_todorov*100 a_todorov*100],'LineWidth',2,'Color',[77 172 38]/255)
xlabel('Mean Force (%Maximum Force)','FontSize',14)
ylabel('CoV (%)','FontSize',14)
xlim([0 100])
legend('New model','Moritz et al. 2005','','Jones et al. 2002','location','northeast')
set(gca,'TickDir','out');
set(gca,'box','off')

figure(2)
shadedErrorBar(mean(mean_Force_norm),mean(std_Force_dt),std(std_Force_dt),'lineprops',{'color',color_code,'LineWidth',2,'markerfacecolor',color_code});
hold on
plot(force_vec,std_vec,'--k','LineWidth',2)
plot(force_vec,std_vec,'o','LineWidth',2,'color','k')
plot(x,x*0.02,'LineWidth',2,'Color',[77 172 38]/255)
plot(x,x*a_todorov,'LineWidth',2,'Color',[208 28 139]/255)
yticks(0:0.25:1.75)
xlabel('Mean Force (%Maximum Force)','FontSize',14)
ylabel('SD (%Maximum Force)','FontSize',14)
ylim([0 1.5])
xlim([0 100])
legend('New model','Moritz et al. 2005','','Jones et al. 2002','Todorov 2004','location','northeast')
set(gca,'TickDir','out');
set(gca,'box','off')
