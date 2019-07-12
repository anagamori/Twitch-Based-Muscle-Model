%==========================================================================
% modelComparison_Tendon.m
% Author: Akira Nagamori
% Last update: 6/27/19
% Descriptions:
%   This code is used to generate fig XX in the paper
%==========================================================================
close all
clear all
clc

%%
condition = 'Model_4_10_CoV_50_Ur_Rec_3';
code_folder = '/Users/akiranagamori/Documents/Github/Twitch-Based-Muscle-Model/Code for Figures';
figure_folder = '/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Figures';

amp_vec = [0.05 0.1:0.1:1];
mean_Force = zeros(10,length(amp_vec));
std_Force = zeros(10,length(amp_vec));
cov_Force = zeros(10,length(amp_vec));
pxx = zeros(10,201);
mean_pxx = zeros(length(amp_vec),201);
f = 0:0.5:100;

for i = 1:4
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
        color_code = [217 4 41]/255;
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
        color_code = [37  65 178]/255;
        vec = [0.1*ones(10,1);0.2*ones(10,1);0.3*ones(10,1);0.4*ones(10,1);0.5*ones(10,1);0.6*ones(10,1);0.7*ones(10,1);0.8*ones(10,1);0.9*ones(10,1);ones(10,1)];
        vec2 = reshape(std_Force,[],1);
        cov_1 = cov_Force;
    elseif i == 3
        Fs = 10000;
        time =0:1/Fs:15;
        data_folder = ['/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Data/withTendon/' condition '_noFV'];
        cd(data_folder)
        load('mean_Force')
        load('std_Force')
        load('cov_Force')
        load('mean_pxx')
        cd(code_folder)
        color_code = [230 97 1]/255;
        
    elseif i == 4
        Fs = 10000;
        time =0:1/Fs:15;
        data_folder = ['/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Data/withTendon/' condition '_shortTendon'];
        cd(data_folder)
        load('mean_Force')
        load('std_Force')
        load('cov_Force')
        load('mean_pxx')
        cd(code_folder)
        color_code = [77 172 38]/255;
        cov_2 = cov_Force;
    end
    mean_mean_Force = mean(mean_Force);
    figure(1)
    %plot([0 amp_vec],[0 mean(mean_Force)]./mean_mean_Force(end),'LineWidth',2,'Color',color_code)
    %plot([0 amp_vec]*100,[0 mean(mean_Force)],'LineWidth',2,'Color',color_code)
    shadedErrorBar([0 amp_vec]*100,[0 mean(mean_Force)],[0 std(mean_Force)],'lineprops',{'color',color_code,'LineWidth',2,'markerfacecolor',color_code})
    hold on
    
    figure(2)
    %errorbar(mean(mean_Force)./mean_mean_Force(end),mean(std_Force),std(std_Force),'LineWidth',2,'Color',color_code);
    %%errorbar(amp_vec*100,mean(std_Force),std(std_Force),'LineWidth',2,'Color',color_code);
    shadedErrorBar(amp_vec*100,mean(std_Force),std(std_Force),'lineprops',{'color',color_code,'LineWidth',2,'markerfacecolor',color_code});
    %shadedErrorBar(mean(mean_Force)/mean_mean_Force(end)*100,mean(std_Force),std(std_Force),'lineprops',{'color',color_code,'LineWidth',2,'markerfacecolor',color_code});
    hold on
    
    
    figure(3)
    %errorbar(mean(mean_Force)./mean_mean_Force(end),mean(cov_Force),std(cov_Force),'LineWidth',2,'Color',color_code);
    %errorbar(amp_vec*100,mean(cov_Force),std(cov_Force),'LineWidth',2,'Color',color_code);
    shadedErrorBar(amp_vec*100,mean(cov_Force),std(cov_Force),'lineprops',{'color',color_code,'LineWidth',2,'markerfacecolor',color_code});
    %shadedErrorBar(mean(mean_Force)/mean_mean_Force(end)*100,mean(cov_Force),std(cov_Force),'lineprops',{'color',color_code,'LineWidth',2,'markerfacecolor',color_code});
    hold on
    
    %%
    figure(4)
    subplot(2,2,1);
    plot(f,mean_pxx(1,:)./sum(mean_pxx(1,:))*100,'LineWidth',2,'color',color_code)
    hold on
    xlim([0 50])
    ylim([0 15])
    yticks(0:5:15)
    xlabel('Frequency (Hz)','FontSize',10)
    ylabel('Power (N^2)','FontSize',10)
    set(gca,'TickDir','out');
    set(gca,'box','off')
    ax = gca;
    ax.FontSize = 6;
    
    
    subplot(2,2,2);
    plot(f,mean_pxx(3,:)./sum(mean_pxx(3,:))*100,'LineWidth',2,'color',color_code)
    hold on
    xlim([0 50])
     ylim([0 15])
    yticks(0:5:15)
    xlabel('Frequency (Hz)','FontSize',10)
    ylabel('Power (N^2)','FontSize',10)
    set(gca,'TickDir','out');
    set(gca,'box','off')
    ax = gca;
    ax.FontSize = 6;
    hold on
    
    subplot(2,2,3);
    plot(f,mean_pxx(5,:)./sum(mean_pxx(5,:))*100,'LineWidth',2,'color',color_code)
    hold on
    xlim([0 50])
    ylim([0 15])
    yticks(0:5:15)
    xlabel('Frequency (Hz)','FontSize',10)
    ylabel('Power (N^2)','FontSize',10)
    set(gca,'TickDir','out');
    set(gca,'box','off')
    ax = gca;
    ax.FontSize = 6;
    
    
    subplot(2,2,4);
    plot(f,mean_pxx(8,:)./sum(mean_pxx(8,:))*100,'LineWidth',2,'color',color_code)
    xlim([0 50])
     ylim([0 15])
    yticks(0:5:15)
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
xlabel('Activation (%)','FontSize',14)
ylabel('Force (N)','FontSize',14)
set(gca,'TickDir','out');
set(gca,'box','off')
legend('Without Tendon','With Tendon','With Tendon & no FV','Shorter Tendon','location','northwest')
% cd (figure_folder)
% saveas(gcf,'activation2meanForce_FV_comparison','pdf')
% cd (code_folder)

figure(2)
%xlabel('Mean Force (%)','FontSize',14)
xlabel('Activation (%)','FontSize',14)
ylabel('SD (N)','FontSize',14)
legend('Without Tendon','With Tendon','With Tendon & no FV','Shorter Tendon','location','northwest')
yticks(0:0.04:0.2)
xlim([0 100])
set(gca,'TickDir','out');
set(gca,'box','off')
% cd (figure_folder)
% saveas(gcf,'activation2SD_FV_comparison','pdf')
% cd (code_folder)

figure(3)
%xlabel('Mean Force (%)','FontSize',14)
xlabel('Activation (%)','FontSize',14)
ylabel('CoV (%)','FontSize',14)
xlim([0 100])
legend('Without Tendon','With Tendon','With Tendon & no FV','Shorter Tendon')
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


