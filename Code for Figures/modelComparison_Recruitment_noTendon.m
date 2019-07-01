%==========================================================================
% analysis_comparison.m
% Author: Akira Nagamori
% Last update: 6/27/19
% Descriptions:
%   This code is used to generate fig. d-f (Summary_Recruitment.pdf)
%==========================================================================
close all
clear all
clc

%%

code_folder = '/Users/akiranagamori/Documents/Github/Twitch-Based-Muscle-Model/Code for Figures';
figure_folder = '/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Figures';

f = 0:0.5:100;

for i = 1:3
    if i == 1
        condition = '10_CoV_50_Ur_Rec_2';
        amp_vec = [0.05 0.1:0.1:1];
        data_folder = ['/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Data/noTendon/' condition];
        cd(data_folder)
        load('mean_Force')
        load('std_Force')
        load('cov_Force')
        load('mean_pxx')
        cd(code_folder)
        color_code = [37  65 178]/255;
    elseif i == 2
        condition = 'Model_4_10_CoV_50_Ur_Rec_3';
        amp_vec = [0.05 0.1:0.1:1];
        data_folder = ['/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Data/noTendon/' condition];
        cd(data_folder)
        load('mean_Force')
        load('std_Force')
        load('cov_Force')
        load('mean_pxx')
        cd(code_folder)
        color_code = [230 57 70]/255;
    elseif i == 3
        condition = '10_CoV_50_Ur_Rec_1';
        amp_vec = [0.05 0.1:0.1:1];
        data_folder = ['/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Data/noTendon/' condition];
        cd(data_folder)
        load('mean_Force')
        load('std_Force')
        load('cov_Force')
        load('mean_pxx')
        cd(code_folder)
        color_code = [100 100 100]/255; 
    end
    mean_mean_Force = mean(mean_Force);
    mean_Force_norm = mean_Force./mean_mean_Force(end)*100;
    figure(1)
    %plot([0 amp_vec],[0 mean(mean_Force)]./mean_mean_Force(end),'LineWidth',2,'Color',color_code)
    %plot([0 amp_vec]*100,[0 mean(mean_Force)],'LineWidth',2,'Color',color_code)
    shadedErrorBar([0 amp_vec]*100,[0 mean(mean_Force_norm)],[0 std(mean_Force_norm)],'lineprops',{'color',color_code,'LineWidth',2,'markerfacecolor',color_code})
    hold on
    
    std_Force_norm = std_Force./mean_mean_Force(end)*100;
    figure(2)
    %shadedErrorBar(amp_vec*100,mean(std_Force_norm),std(std_Force_norm),'lineprops',{'color',color_code,'LineWidth',2,'markerfacecolor',color_code});
    shadedErrorBar(mean(mean_Force_norm),mean(std_Force_norm),std(std_Force_norm),'lineprops',{'color',color_code,'LineWidth',2,'markerfacecolor',color_code});
    hold on
    
    
    figure(3)
    %shadedErrorBar(amp_vec*100,mean(cov_Force),std(cov_Force),'lineprops',{'color',color_code,'LineWidth',2,'markerfacecolor',color_code});
    shadedErrorBar(mean(mean_Force_norm),mean(cov_Force),std(cov_Force),'lineprops',{'color',color_code,'LineWidth',2,'markerfacecolor',color_code});
    hold on
    
    %%
    figure(4)
    subplot(2,2,1);
    plot(f,mean_pxx(1,:),'LineWidth',2,'color',color_code)
    hold on
    xlim([0 50])
    xlabel('Frequency (Hz)','FontSize',10)
    ylabel('Power (N^2)','FontSize',10)
    set(gca,'TickDir','out');
    set(gca,'box','off')
    ax = gca;
    ax.FontSize = 6;
    
    
    subplot(2,2,2);
    plot(f,mean_pxx(3,:),'LineWidth',2,'color',color_code)
    hold on
    xlim([0 50])
    xlabel('Frequency (Hz)','FontSize',10)
    ylabel('Power (N^2)','FontSize',10)
    set(gca,'TickDir','out');
    set(gca,'box','off')
    ax = gca;
    ax.FontSize = 6;
    hold on
    
    subplot(2,2,3);
    plot(f,mean_pxx(5,:),'LineWidth',2,'color',color_code)
    hold on
    xlim([0 50])
    xlabel('Frequency (Hz)','FontSize',10)
    ylabel('Power (N^2)','FontSize',10)
    set(gca,'TickDir','out');
    set(gca,'box','off')
    ax = gca;
    ax.FontSize = 6;
    
    
    subplot(2,2,4);
    plot(f,mean_pxx(8,:),'LineWidth',2,'color',color_code)
    xlim([0 50])
    hold on
    xlabel('Frequency (Hz)','FontSize',10)
    ylabel('Power (N^2)','FontSize',10)
    set(gca,'TickDir','out');
    set(gca,'box','off')
    ax = gca;
    ax.FontSize = 6;
    
end

%%
figure(1)
xlabel('Activation (%Maxmum)','FontSize',14)
ylabel('Force (%Maximum)','FontSize',14)
ylim([0 105])
set(gca,'TickDir','out');
set(gca,'box','off')
legend('Constant','Combined','Variable','location','northwest')
% cd (figure_folder)
% saveas(gcf,'activation2meanForce_FV_comparison','pdf')
% cd (code_folder)

figure(2)
%xlabel('Mean Force (%)','FontSize',14)
xlabel('Mean Force (%Maximum Force)','FontSize',14)
ylabel('SD (%Maximum Force)','FontSize',14)
legend('Constant','Combined','Variable','location','northwest')
%yticks([0.05 0.1 0.15 0.2 0.25])
xlim([0 100])
set(gca,'TickDir','out');
set(gca,'box','off')
% cd (figure_folder)
% saveas(gcf,'activation2SD_FV_comparison','pdf')
% cd (code_folder)

figure(3)
%xlabel('Mean Force (%)','FontSize',14)
xlabel('Mean Force (%Maximum Force)','FontSize',14)
ylabel('CoV (%)','FontSize',14)
xlim([0 100])
legend('Constant','Combined','Variable','location','northeast')
set(gca,'TickDir','out');
set(gca,'box','off')
% cd (figure_folder)
% saveas(gcf,'activation2CoV_FV_comparison','pdf')
% cd (code_folder)

figure(4)
fig = gcf;
%linkaxes([ax1,ax2,ax3,ax4],'y')
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 4.56 4.56];
% cd (figure_folder)
% saveas(gcf,'pxx_FV_comparison','pdf')
% cd (code_folder)

