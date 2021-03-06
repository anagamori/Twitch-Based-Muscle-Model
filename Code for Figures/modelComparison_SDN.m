%==========================================================================
% modelComparison_SDN.m
% Author: akiranagamori Nagamori
% Last update: 6/27/19
% Descriptions:
%   This code is used to generate results_SDN.pdf
%==========================================================================
close all
clear all
clc

%%
code_folder = '/Users/akira/Documents/Github/Twitch-Based-Muscle-Model/Code for Figures';
figure_folder = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Figures';

amp_vec = [0.025 0.05 0.1:0.1:1];
mean_Force = zeros(10,length(amp_vec));
std_Force = zeros(10,length(amp_vec));
cov_Force = zeros(10,length(amp_vec));
pxx = zeros(10,201);
mean_pxx = zeros(length(amp_vec),201);
f = 0:0.5:100;

for i = 1:8
    if i == 1
        condition = 'Model_default_v2'; %_CTvsPTi';
        Fs = 2000;
        time =0:1/Fs:15;
        data_folder = ['/Users/akira/Documents/Github/Twitch-Based-Muscle-Model/Data/New Model//' condition];
        cd(data_folder)
        load('mean_Force')
        load('std_Force')
        load('std_Force_dt')
        load('cov_Force')
        load('cov_Force_dt')
        load('mean_pxx')
        cd(code_folder)
        color_code = [0 0 0]/255;
        mean_mean_Force = mean(mean_Force);
        MVC_default = mean_mean_Force(end);
        std_default = mean(std_Force./MVC_default.*100);
        mean_force_default = mean(mean_Force)./mean_mean_Force(end)*100;
    elseif i == 2
        condition = 'Model_PR_100'; %_CTvsPTi_PR_100';
        Fs = 2000;
        time =0:1/Fs:15;
        data_folder = ['/Users/akira/Documents/Github/Twitch-Based-Muscle-Model/Data/New Model/' condition];
        cd(data_folder)
        load('mean_Force')
        load('std_Force')
        load('std_Force_dt')
        load('cov_Force')
        load('cov_Force_dt')
        load('mean_pxx')
        cd(code_folder)
        color_code = [255 186 8]/255;
        %color_code = [77 172 38]/255;
    elseif i == 3
        condition = 'Model_Ur_50';
        Fs = 2000;
        time =0:1/Fs:15;
        data_folder = ['/Users/akira/Documents/Github/Twitch-Based-Muscle-Model/Data/New Model/' condition];
        cd(data_folder)
        load('mean_Force')
        load('std_Force')
        load('std_Force_dt')
        load('cov_Force')
        load('cov_Force_dt')
        load('mean_pxx')
        cd(code_folder)
        %color_code =  [208 28 139]/255;
        color_code = [19 111 99]/255;
        vec = [0.1*ones(10,1);0.2*ones(10,1);0.3*ones(10,1);0.4*ones(10,1);0.5*ones(10,1);0.6*ones(10,1);0.7*ones(10,1);0.8*ones(10,1);0.9*ones(10,1);ones(10,1)];
        vec2 = reshape(std_Force,[],1);
         elseif i == 4
        condition = 'Model_N_100';
        Fs = 2000;
        time =0:1/Fs:15;
        data_folder = ['/Users/akira/Documents/Github/Twitch-Based-Muscle-Model/Data/New Model/' condition];
        cd(data_folder)
        load('mean_Force')
        load('std_Force')
        load('std_Force_dt')
        load('cov_Force')
        load('cov_Force_dt')
        load('mean_pxx')
        cd(code_folder)
        color_code = [76 201 240]/255;
        
        vec = [0.1*ones(10,1);0.2*ones(10,1);0.3*ones(10,1);0.4*ones(10,1);0.5*ones(10,1);0.6*ones(10,1);0.7*ones(10,1);0.8*ones(10,1);0.9*ones(10,1);ones(10,1)];
        vec2 = reshape(std_Force,[],1);
    elseif i == 5
        condition = 'Model_onion_skin_v4';
        Fs = 2000;
        time =0:1/Fs:15;
        data_folder = ['/Users/akira/Documents/Github/Twitch-Based-Muscle-Model/Data/New Model/' condition];
        cd(data_folder)
        load('mean_Force')
        load('std_Force')
        load('std_Force_dt')
        load('cov_Force')
        load('cov_Force_dt')
        load('mean_pxx')
        cd(code_folder)
        color_code = [100 100 100]/255;
        vec = [0.1*ones(10,1);0.2*ones(10,1);0.3*ones(10,1);0.4*ones(10,1);0.5*ones(10,1);0.6*ones(10,1);0.7*ones(10,1);0.8*ones(10,1);0.9*ones(10,1);ones(10,1)];
        vec2 = reshape(std_Force,[],1);
    elseif i == 6
        condition = 'Model_default_v2';
        Fs = 2000;
        time =0:1/Fs:15;
        data_folder = ['/Users/akira/Documents/Github/Twitch-Based-Muscle-Model/Data/New Model/noTendon/' condition];
        cd(data_folder)
        load('mean_Force')
        load('std_Force')
        load('std_Force_dt')
        load('cov_Force')
        load('cov_Force_dt')
        load('mean_pxx')
        cd(code_folder)
        color_code = [230 57 70]/255;
        vec = [0.1*ones(10,1);0.2*ones(10,1);0.3*ones(10,1);0.4*ones(10,1);0.5*ones(10,1);0.6*ones(10,1);0.7*ones(10,1);0.8*ones(10,1);0.9*ones(10,1);ones(10,1)];
        vec2 = reshape(std_Force,[],1);
     elseif i == 7
        condition = 'Model_onion_skin_N_100_Ur_80';
        Fs = 2000;
        time =0:1/Fs:15;
        data_folder = ['/Users/akira/Documents/Github/Twitch-Based-Muscle-Model/Data/New Model/' condition];
        cd(data_folder)
        load('mean_Force')
        load('std_Force')
        load('std_Force_dt')
        load('cov_Force')
        load('cov_Force_dt')
        load('mean_pxx')
        cd(code_folder)
        color_code = [114 9 183]/255;
        vec = [0.1*ones(10,1);0.2*ones(10,1);0.3*ones(10,1);0.4*ones(10,1);0.5*ones(10,1);0.6*ones(10,1);0.7*ones(10,1);0.8*ones(10,1);0.9*ones(10,1);ones(10,1)];
        vec2 = reshape(std_Force,[],1);    
     elseif i == 8
        condition = 'Model_onion_skin_N_100_Ur_80';
        Fs = 2000;
        time =0:1/Fs:15;
        data_folder = ['/Users/akira/Documents/Github/Twitch-Based-Muscle-Model/Data/New Model/noTendon/' condition];
        cd(data_folder)
        load('mean_Force')
        load('std_Force')
        load('std_Force_dt')
        load('cov_Force')
        load('cov_Force_dt')
        load('mean_pxx')
        cd(code_folder)
        %color_code = [76 201 240]/255;
        color_code = [230 97 1]/255;
        vec = [0.1*ones(10,1);0.2*ones(10,1);0.3*ones(10,1);0.4*ones(10,1);0.5*ones(10,1);0.6*ones(10,1);0.7*ones(10,1);0.8*ones(10,1);0.9*ones(10,1);ones(10,1)];
        vec2 = reshape(std_Force,[],1);    
    
    end
    mean_mean_Force = mean(mean_Force);
    MVC = mean_mean_Force(end);
    figure(1)
    %plot([0 amp_vec],[0 mean(mean_Force)]./mean_mean_Force(end),'LineWidth',2,'Color',color_code)
    %plot([0 amp_vec]*100,[0 mean(mean_Force)],'LineWidth',2,'Color',color_code)
    %shadedErrorBar([0 amp_vec]*100,[0 mean(mean_Force./MVC.*100)],[0 std(mean_Force./MVC.*100)],'lineprops',{'color',color_code,'LineWidth',2,'markerfacecolor',color_code})
    shadedErrorBar([0 amp_vec]*100,[0 mean(mean_Force)],[0 std(mean_Force)],'lineprops',{'color',color_code,'LineWidth',2,'markerfacecolor',color_code})
    hold on
    
    figure(2)
    %errorbar(mean(mean_Force)./mean_mean_Force(end),mean(std_Force),std(std_Force),'LineWidth',2,'Color',color_code);
    %errorbar(amp_vec*100,mean(std_Force),std(std_Force),'LineWidth',2,'Color',color_code);
    shadedErrorBar(mean(mean_Force)./mean_mean_Force(end)*100,mean(std_Force_dt./MVC.*100),std(std_Force_dt./MVC.*100),'lineprops',{'color',color_code,'LineWidth',2,'markerfacecolor',color_code});
    %shadedErrorBar(amp_vec*100,mean(std_Force./MVC.*100),std(std_Force./MVC.*100),'lineprops',{'color',color_code,'LineWidth',2,'markerfacecolor',color_code});
    hold on
    
    
    figure(3)
    %errorbar(mean(mean_Force)./mean_mean_Force(end),mean(cov_Force),std(cov_Force),'LineWidth',2,'Color',color_code);
    %errorbar(amp_vec*100,mean(cov_Force),std(cov_Force),'LineWidth',2,'Color',color_code);
    %shadedErrorBar(amp_vec*100,mean(cov_Force),std(cov_Force),'lineprops',{'color',color_code,'LineWidth',2,'markerfacecolor',color_code});
    shadedErrorBar(mean(mean_Force)./mean_mean_Force(end)*100,mean(cov_Force_dt),std(cov_Force_dt),'lineprops',{'color',color_code,'LineWidth',2,'markerfacecolor',color_code});
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
xlabel('Activation (%)','FontSize',14)
ylabel('Force (% Maximum)','FontSize',14)
set(gca,'TickDir','out');
set(gca,'box','off')
%legend('Default','RP = 100','Ur = 0.5','N = 400','location','northwest')
%ylim([0 101])
% cd (figure_folder)
% saveas(gcf,'activation2meanForce_FV_comparison','pdf')
% cd (code_folder)

a = std_default(1)/mean_force_default(1);
a_todorov = 0.1289;
x = [0 amp_vec]*100;
% y_0 = 0.2*1/6;
% y_100 = 1.2+0.2*2.5/6;
% a = (y_100-y_0)/100;
% 
% std_N_100 = mean(std_Force./MVC.*100);
% mean_force_N_100 = mean(mean_Force)./mean_mean_Force(end)*100;
% y_5_new = std_N_100(1);
% b = y_5_new - a*mean_force_N_100(1);

figure(2)
xlabel('Mean Force (%)','FontSize',14)
plot(x,0.0186*x.^1.05,'--','LineWidth',1,'Color','k')
plot(x,x*a_todorov,'LineWidth',1,'Color','k')
hold on 
%plot(x,x*a_todorov,'--','LineWidth',2,'Color','k')
xlabel('Mean Force (%)','FontSize',14)
ylabel('SD (%MVC)','FontSize',14)
legend('Default','RP = 100','Ur = 0.5','N = 100','Oion-Skin','no SEE','Onion-skin + SSE + RP = 100','Onion-skin + no SSE + RP = 100','Jones et al. 2002','location','northwest')
yticks(0:0.5:2.5)
ylim([0 2.5])
xlim([0 100])
set(gca,'TickDir','out');
set(gca,'box','off')
% cd (figure_folder)
% saveas(gcf,'activation2SD_FV_comparison','pdf')
% cd (code_folder)

y = 0.0186*x.^1.05./x*100;
figure(3)
%xlabel('Mean Force (%)','FontSize',14)
plot(x,y,'--k','LineWidth',1)
xlabel('Activation (%)','FontSize',14)
ylabel('CoV (%)','FontSize',14)
xlim([0 100])
legend('Default','RP = 100','Ur = 0.5','N = 400')
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
% saveas(gcf,'pxx_SDN_comparison','pdf')
% cd (code_folder)

