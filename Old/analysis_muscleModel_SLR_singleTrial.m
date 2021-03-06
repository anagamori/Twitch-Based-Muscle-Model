%==========================================================================
% analysis_muscleModel.m
% Author: Akira Nagamori
% Last update: 3/5/19
% Descriptions:
%   Run muscle model simulation
%==========================================================================
close all
clear all
clc

%%
condition = 'GD_40_GS_40_Ia_2000_DL_30';
data_folder = ['/Volumes/DATA2/New_Model/SLR/' condition];
code_folder = '/Users/akira/Documents/Github/Twitch-Based-Muscle-Model';
figure_folder = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Figures';

%%
amp_vec = [0.05 0.1:0.1:1];

% %%
for j = 1  %:length(amp_vec)
    if j < 2
        Fs = 10000;
        time = 0:1/Fs:15;
    elseif j >= 2 && j <= 3
        Fs = 15000;
        time = 0:1/Fs:15;
    elseif j >= 4 && j < 7
        Fs = 20000;
        time = 0:1/Fs:15;
    elseif j >= 7
        Fs = 25000;
        time = 0:1/Fs:15;
    end
%         if j <= 1
%             Fs = 10000;
%             time = 0:1/Fs:15;
%         elseif j >= 2 && j <= 3
%             Fs = 15000;
%             time = 0:1/Fs:15;
%         elseif j == 4
%             Fs = 20000;
%             time = 0:1/Fs:15;
%         elseif j >= 5 && j < 7
%             Fs = 25000;
%             time = 0:1/Fs:15;
%         elseif j >= 7
%             Fs = 30000;
%             time = 0:1/Fs:15;
%         end
%     
    j
    tic
    for i = 7
        cd(data_folder)
        load(['Data_' num2str(j) '_' num2str(i)])
        cd(code_folder)
        Force =  output.ForceTendon;
        mean_Force(j+1) = mean(Force(5*Fs+1:end));
        std_Force(j+1) = std(Force(5*Fs+1:end));
        cov_Force(j+1) =  std_Force(j+1)/mean_Force(j+1)*100;
        [pxx,f] = pwelch(Force(5*Fs+1:end)-mean(Force(5*Fs+1:end)),gausswin(5*Fs),0.9*5*Fs,0:0.1:100,Fs,'power');
        
        figure(11)
        plot(time,Force)
        hold on
        
        figure(12)
        plot(f,pxx)
        xlim([0 50])
        
    end
    toc
    mean_pxx(j+1,:) = pxx;

end

%%

mean_mean_Force = mean(mean_Force);
figure(1)
plot([0 amp_vec(1:j+1)],[0 mean_Force],'LineWidth',2)
%plot([0 amp_vec],[0 mean_Force]./mean_Force(end),'LineWidth',2)
xlabel('Activation (%)','FontSize',14)
ylabel('Force (%MVC)','FontSize',14)
set(gca,'TickDir','out');
set(gca,'box','off')
% cd (figure_folder)
% saveas(gcf,'activation2meanForce_withTendon','pdf')
% cd (code_folder)

figure(2)
plot(amp_vec(1:j+1),std_Force,'LineWidth',2)
%plot(amp_vec,std_Force,'LineWidth',2)
xlabel('Mean Force (%)','FontSize',14)
ylabel('SD (N)','FontSize',14)
set(gca,'TickDir','out');
set(gca,'box','off')
% cd (figure_folder)
% saveas(gcf,'meanForce2SD_withTendon','pdf')
% cd (code_folder)

figure(3)
%plot(amp_vec,cov_Force,'LineWidth',2)
plot(amp_vec(1:j+1),cov_Force,'LineWidth',2)
% hold on
% plot(mean(mean_Force)./mean_mean_Force(end),mean(cov_Force),'LineWidth',2)
xlabel('Mean Force (%)','FontSize',14)
ylabel('CoV (%)','FontSize',14)
set(gca,'TickDir','out');
set(gca,'box','off')
% cd (figure_folder)
% saveas(gcf,'meanForce2CoV_withTendon','pdf')
% cd (code_folder)

figure(4)
plot(f,mean_pxx,'LineWidth',2)
xlabel('Frequency (Hz)','FontSize',14)
ylabel('Power (N^2)','FontSize',14)
set(gca,'TickDir','out');
set(gca,'box','off')
legend('10% MVC','30% MVC','50% MVC','80% MVC')
% cd (figure_folder)
% saveas(gcf,'pxx_withTendon','pdf')
% cd (code_folder)