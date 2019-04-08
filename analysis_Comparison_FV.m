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

for i = 1:3
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
        color_code = [77 172 38]/255;
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
        color_code = [123 50 148]/255;
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
%     elseif i == 4
%         Fs = 10000;
%         time =0:1/Fs:15;
%         data_folder = ['/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Data/withTendon/' condition '_shortTendon'];
%         cd(data_folder)
%         load('mean_Force')
%         load('std_Force')
%         load('cov_Force')
%         load('mean_pxx')
%         cd(code_folder)
%         color_code = [208 28 139]/255;
    end
    mean_mean_Force = mean(mean_Force);
    figure(1)
    plot([0 amp_vec],[0 mean(mean_Force)]./mean_mean_Force(end),'LineWidth',2,'Color',color_code)
    %plot([0 amp_vec],[0 mean(mean_Force)],'LineWidth',2)
    xlabel('Activation (%)','FontSize',14)
    ylabel('Force (%MVC)','FontSize',14)
    set(gca,'TickDir','out');
    set(gca,'box','off')
    hold on
    
    figure(2)
    errorbar(mean(mean_Force)./mean_mean_Force(end),mean(std_Force),std(std_Force),'LineWidth',2,'Color',color_code);
    %errorbar(amp_vec,mean(std_Force),std(std_Force),'LineWidth',2);
    xlabel('Mean Force (%)','FontSize',14)
    ylabel('SD (N)','FontSize',14)
    hold on
    yticks([0.05 0.1 0.15 0.2 0.25])
    set(gca,'TickDir','out');
    set(gca,'box','off')

    figure(3)
    errorbar(mean(mean_Force)./mean_mean_Force(end),mean(cov_Force),std(cov_Force),'LineWidth',2,'Color',color_code);
    %errorbar(amp_vec,mean(cov_Force),std(cov_Force),'LineWidth',2);
    xlabel('Mean Force (%)','FontSize',14)
    ylabel('CoV (%)','FontSize',14)
    hold on
    set(gca,'TickDir','out');
    set(gca,'box','off')
    
    %%
    figure(4)
    if i == 1
        ax1 = subplot(1,3,1);
    elseif i == 2
        ax2 = subplot(1,3,2);
    elseif i == 3
        ax3 = subplot(1,3,3);
    end
    plot(f,mean_pxx([1 3 5 8],:),'LineWidth',2)
    xlim([0 50])
    xlabel('Frequency (Hz)','FontSize',14)
    ylabel('Power (N^2)','FontSize',14)
    set(gca,'TickDir','out');
    set(gca,'box','off')
    legend('10% MVC','30% MVC','50% MVC','80% MVC')
    
end

figure(1)
legend('Without Tendon','With Tendon','With Tendon & no FV','location','northwest')
cd (figure_folder)
saveas(gcf,'activation2meanForce_FV_comparison','pdf')
cd (code_folder)

figure(2)
legend('Without Tendon','With Tendon','With Tendon & no FV','location','northwest')
cd (figure_folder)
saveas(gcf,'meanForce2SD_FV_comparison','pdf')
cd (code_folder)

figure(3)
legend('Without Tendon','With Tendon','With Tendon & no FV')
cd (figure_folder)
saveas(gcf,'meanForce2CoV_FV_comparison','pdf')
cd (code_folder)

figure(4)
linkaxes([ax1,ax2,ax3],'y')
cd (figure_folder)
saveas(gcf,'pxx_FV_comparison','pdf')
cd (code_folder)

