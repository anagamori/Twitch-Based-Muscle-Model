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
condition = 'Model_4_20_CoV_50_Ur_Rec_3';
code_folder = '/Users/akira/Documents/Github/Twitch-Based-Muscle-Model/Code for Figures';
figure_folder = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Figures';

f = 0:0.5:100;

for i = 1:3
    if i == 1
        amp_vec = 0.1:0.1:1;
        data_folder = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Data/Fuglevand/Model_1_CoV_20';
        cd(data_folder)
        load('mean_Force')
        load('std_Force')
        load('cov_Force')
        load('mean_pxx')
        cd(code_folder)
        color_code = [100 100 100]/255;
    elseif i == 2
        amp_vec = [0.05 0.1:0.1:1];
        data_folder = ['/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Data/noTendon/' condition];
        cd(data_folder)
        load('mean_Force')
        load('std_Force')
        load('cov_Force')
        load('mean_pxx')
        cd(code_folder)
        color_code = [230 57 70]/255;
    elseif i == 3
        amp_vec = [0.05 0.1:0.1:1];
        data_folder = ['/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Data/withTendon/' condition];
        cd(data_folder)
        load('mean_Force')
        load('std_Force')
        load('cov_Force')
        load('mean_pxx')
        cd(code_folder)
        color_code = [37  65 178]/255;
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
    plot(f,mean_pxx(1,:)./sum(mean_pxx(1,:))*100,'LineWidth',2,'color',color_code)
    hold on
    xlim([0 50])
    ylim([0 20])
    yticks(0:5:25)
    xlabel('Frequency (Hz)','FontSize',10)
    ylabel('Proportion of Total Power (%)','FontSize',10)
    set(gca,'TickDir','out');
    set(gca,'box','off')
    ax = gca;
    ax.FontSize = 6;
    
    
    subplot(2,2,2);
    plot(f,mean_pxx(3,:)./sum(mean_pxx(3,:))*100,'LineWidth',2,'color',color_code)
    hold on
    xlim([0 50])
    ylim([0 20])
    yticks(0:5:25)
    xlabel('Frequency (Hz)','FontSize',10)
    ylabel('Proportion of Total Power (%)','FontSize',10)
    set(gca,'TickDir','out');
    set(gca,'box','off')
    ax = gca;
    ax.FontSize = 6;
    hold on
    
    subplot(2,2,3);
    plot(f,mean_pxx(5,:)./sum(mean_pxx(5,:))*100,'LineWidth',2,'color',color_code)
    hold on
    xlim([0 50])
    ylim([0 20])
    yticks(0:5:25)
    xlabel('Frequency (Hz)','FontSize',10)
    ylabel('Proportion of Total Power (%)','FontSize',10)
    set(gca,'TickDir','out');
    set(gca,'box','off')
    ax = gca;
    ax.FontSize = 6;
    
    
    subplot(2,2,4);
    plot(f,mean_pxx(8,:)./sum(mean_pxx(8,:))*100,'LineWidth',2,'color',color_code)
    xlim([0 50])
    ylim([0 20])
    yticks(0:5:25)
    hold on
    xlabel('Frequency (Hz)','FontSize',10)
    ylabel('Proportion of Total Power (%)','FontSize',10)
    set(gca,'TickDir','out');
    set(gca,'box','off')
    ax = gca;
    ax.FontSize = 6;
    
end

%%
figure(1)
xlabel('Activation (%)','FontSize',14)
ylabel('Force (%Maximum)','FontSize',14)
ylim([0 105])
set(gca,'TickDir','out');
set(gca,'box','off')
legend('Fuglevand model','New model without tendon','New model with tendon','location','northwest')
% cd (figure_folder)
% saveas(gcf,'activation2meanForce_FV_comparison','pdf')
% cd (code_folder)

figure(2)
%xlabel('Mean Force (%)','FontSize',14)
xlabel('Mean Force (%Maximum Force)','FontSize',14)
ylabel('SD (%Maximum Force)','FontSize',14)
legend('Fuglevand model','New model without tendon','New model with tendon','location','northwest')
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
legend('Fuglevand model','New model without tendon','New model without tendon','location','northwest')
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

