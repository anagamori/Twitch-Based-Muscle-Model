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
condition = 'Model_default_v2';
code_folder = '/Users/akiranagamori/Documents/Github/Twitch-Based-Muscle-Model/Code for Figures';
figure_folder = '/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Figures';

f = 0:0.5:100;

for i = 1:2
    if i == 1
        amp_vec = [0.025 0.05  0.1:0.1:1];
        data_folder = '/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Data/Fuglevand/N_200_CoV_var';
        cd(data_folder)
        load('mean_Force')
        load('std_Force')
        load('cov_Force')
        load('cov_Force_dt')
        load('mean_pxx')
        cd(code_folder)
        color_code = [100 100 100]/255;
        mean_mean_Force = mean(mean_Force);
        MVC_FM = mean_mean_Force(end);
        std_FM = mean(std_Force./MVC_FM.*100);
        mean_force_FM = mean(mean_Force)./mean_mean_Force(end)*100;

    elseif i == 2
        amp_vec = [0.025 0.05 0.1:0.1:1];
        data_folder = ['/Users/akiranagamori/Documents/Github/Twitch-Based-Muscle-Model/Data/New Model/' condition];
        cd(data_folder)
        load('mean_Force')
        load('std_Force')
        load('cov_Force')
        load('cov_Force_dt')
        load('mean_pxx')
        cd(code_folder)
        color_code = [230 57 70]/255;
    elseif i == 3
        amp_vec = [0.05 0.1:0.1:1];
        data_folder = ['/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Data/withTendon/' condition];
        cd(data_folder)
        load('mean_Force')
        load('std_Force')
        load('cov_Force')
        load('cov_Force_dt')
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
    shadedErrorBar(amp_vec*100,mean(std_Force_norm),std(std_Force_norm),'lineprops',{'color',color_code,'LineWidth',2,'markerfacecolor',color_code});
    %shadedErrorBar(mean(mean_Force_norm),mean(std_Force_norm),std(std_Force_norm),'lineprops',{'color',color_code,'LineWidth',2,'markerfacecolor',color_code});
    hold on
    
    
    figure(3)
    shadedErrorBar(amp_vec*100,mean(cov_Force),std(cov_Force),'lineprops',{'color',color_code,'LineWidth',2,'markerfacecolor',color_code});
    %shadedErrorBar(mean(mean_Force_norm),mean(cov_Force_dt),std(cov_Force_dt),'lineprops',{'color',color_code,'LineWidth',2,'markerfacecolor',color_code});
    hold on
    
     figure(4)
    %shadedErrorBar(amp_vec*100,mean(cov_Force_dt),std(cov_Force_dt),'lineprops',{'color',color_code,'LineWidth',2,'markerfacecolor',color_code});
    shadedErrorBar(mean(mean_Force_norm),mean(cov_Force_dt),std(cov_Force_dt),'lineprops',{'color',color_code,'LineWidth',2,'markerfacecolor',color_code});
    hold on
    %%
    figure(5)
    subplot(1,4,1);
    plot(f,mean_pxx(1,:)./sum(mean_pxx(1,:))*100,'LineWidth',2,'color',color_code)
    hold on
    xlim([0 30])
    ylim([0 20])
    yticks(0:5:25)
    xlabel('Frequency (Hz)','FontSize',10)
    ylabel('Proportion of Total Power (%)','FontSize',10)
    set(gca,'TickDir','out');
    set(gca,'box','off')
    ax = gca;
    ax.FontSize = 6;
    
    
    subplot(1,4,2);
    plot(f,mean_pxx(3,:)./sum(mean_pxx(3,:))*100,'LineWidth',2,'color',color_code)
    hold on
    xlim([0 30])
    ylim([0 20])
    yticks(0:5:25)
    xlabel('Frequency (Hz)','FontSize',10)
    %ylabel('Proportion of Total Power (%)','FontSize',10)
    set(gca,'TickDir','out');
    set(gca,'box','off')
    ax = gca;
    ax.FontSize = 6;
    hold on
    
    subplot(1,4,3);
    plot(f,mean_pxx(5,:)./sum(mean_pxx(5,:))*100,'LineWidth',2,'color',color_code)
    hold on
    xlim([0 30])
    ylim([0 20])
    yticks(0:5:25)
    xlabel('Frequency (Hz)','FontSize',10)
    %ylabel('Proportion of Total Power (%)','FontSize',10)
    set(gca,'TickDir','out');
    set(gca,'box','off')
    ax = gca;
    ax.FontSize = 6;
    
    
    subplot(1,4,4);
    plot(f,mean_pxx(8,:)./sum(mean_pxx(8,:))*100,'LineWidth',2,'color',color_code)
    xlim([0 30])
    ylim([0 20])
    yticks(0:5:25)
    hold on
    xlabel('Frequency (Hz)','FontSize',10)
    %ylabel('Proportion of Total Power (%)','FontSize',10)
    set(gca,'TickDir','out');
    set(gca,'box','off')
    ax = gca;
    ax.FontSize = 6;
    
end

%%
x = [0 amp_vec]*100;
y_0 = 0.2*1/6;
y_100 = 1.2+0.2*2.5/6;
a = (y_100-y_0)/100;

y_5_new = std_FM(1);
b = y_5_new - a*mean_force_FM(1);

force_vec = [2.5 5 10 30 50 80];
CoV_vec = [3.9402 2.7496 1.6535 1.2189 1.4929 2.211];
std_vec = CoV_vec./100.*force_vec;
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
%plot(x,x*a+b,'LineWidth',2,'Color','k')
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
plot(force_vec,CoV_vec,'--k','LineWidth',2)
plot(force_vec,CoV_vec,'o','LineWidth',2,'color','k')

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
%xlabel('Mean Force (%)','FontSize',14)
plot(force_vec,CoV_vec,'--k','LineWidth',2)
plot(force_vec,CoV_vec,'o','LineWidth',2,'color','k')

xlabel('Mean Force (%Maximum Force)','FontSize',14)
ylabel('CoV (%)','FontSize',14)
xlim([0 100])
legend('Fuglevand model','New model without tendon','New model without tendon','location','northwest')
set(gca,'TickDir','out');
set(gca,'box','off')


figure(5)
fig = gcf;
%linkaxes([ax1,ax2,ax3,ax4],'y')
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6.29 4];
print('-fillpage','FillPageFigure','-dpdf')
% cd (figure_folder)
% saveas(gcf,'pxx_FV_comparison','pdf')
% cd (code_folder)

