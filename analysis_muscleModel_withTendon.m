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
condition = 'Model_9_var_CoV_50_Ur_Rec_3';
data_folder = ['/Volumes/DATA2/New_Model/withTendon/' condition];
data_folder_git = ['/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Data/withTendon/' condition];
code_folder = '/Users/akira/Documents/Github/Twitch-Based-Muscle-Model';
figure_folder = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Figures';

%%
amp_vec = [0.05 0.1:0.1:1];

nTrial = 10;
Force_mat = zeros(nTrial,15*10000+1);
mean_Force = zeros(nTrial,length(amp_vec));
std_Force = zeros(nTrial,length(amp_vec));
cov_Force = zeros(nTrial,length(amp_vec));
pxx = zeros(nTrial,1001);
mean_pxx = zeros(length(amp_vec),1001);
mean_pxx_2 = zeros(length(amp_vec),1001);
% %%
for j = 6:7 %:length(amp_vec)
    if j <= 1
        Fs = 10000;
        time = 0:1/Fs:15;
    elseif j >= 2 && j <= 4
        Fs = 15000;
        time = 0:1/Fs:15;
    elseif j > 4 && j < 7
        Fs = 20000;
        time = 0:1/Fs:15;
    elseif j >= 7 && j <= 8
        Fs = 25000;
        time = 0:1/Fs:15;
    elseif j >= 9
        Fs = 30000;
        time = 0:1/Fs:15;
    end
    %             time = 0:1/Fs:15;
    %         elseif j >= 2 && j < 3
    %             Fs = 15000;
    %             time = 0:1/Fs:15;
    %         elseif j >= 3 && j <= 4
    %             Fs = 20000;
    %             time = 0:1/Fs:15;
    %         elseif j >= 5 && j < 6
    %             Fs = 25000;
    %             time = 0:1/Fs:15;
    %         elseif j >= 6 && j<7
    %             Fs = 30000;
    %             time = 0:1/Fs:15;
    %         else
    %             Fs = 40000;
    %             time = 0:1/Fs:15;
    %         end
    %
    j
    tic
    for i = 1:10
        cd(data_folder)
        load(['Data_' num2str(j) '_' num2str(i)])
        cd(code_folder)
        Force =  output.ForceTendon;
        mean_Force(i,j+2) = mean(Force(5*Fs+2:end));
        std_Force(i,j+2) = std(Force(5*Fs+2:end));
        cov_Force(i,j+2) =  std_Force(i,j+2)/mean_Force(i,j+2)*100;
        [pxx(i,:),f] = pwelch(Force(5*Fs+1:end)-mean(Force(5*Fs+1:end)),gausswin(5*Fs),0.9*5*Fs,0:0.1:100,Fs,'power');
        figure(11)
        plot(time,Force)
        hold on
        
        
        idx = linspace(1,length(Force),10000*15+1);
        if j < 2
            Force_mat(i,:) = Force;
        else
            Force_mat(i,:) = interp1(1:length(Force),Force,idx,'linear');
        end
        %Force_mat(i,:) = Force;
    end
    toc
    mean_pxx(j+2,:) = mean(pxx);
    cd(data_folder)
    save(['Force_mat_' num2str(j)],'Force_mat')
    cd(code_folder)
    cd(data_folder_git)
    save(['Force_mat_' num2str(j)],'Force_mat')
    cd(code_folder)
end

%%

mean_mean_Force = mean(mean_Force);
figure(1)
plot([0 amp_vec],[0 mean(mean_Force)]./mean_mean_Force(end),'LineWidth',2)
%plot([0 amp_vec],[0 mean_Force]./mean_Force(end),'LineWidth',2)
xlabel('Activation (%)','FontSize',14)
ylabel('Force (%MVC)','FontSize',14)
set(gca,'TickDir','out');
set(gca,'box','off')
% cd (figure_folder)
% saveas(gcf,'activation2meanForce_withTendon','pdf')
% cd (code_folder)

figure(2)
plot(amp_vec,mean(std_Force),'LineWidth',2)
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
plot(amp_vec,mean(cov_Force),'LineWidth',2)
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
plot(f,mean_pxx([1 3 5 8],:),'LineWidth',2)
xlabel('Frequency (Hz)','FontSize',14)
ylabel('Power (N^2)','FontSize',14)
set(gca,'TickDir','out');
set(gca,'box','off')
legend('10% MVC','30% MVC','50% MVC','80% MVC')
% cd (figure_folder)
% saveas(gcf,'pxx_withTendon','pdf')
% cd (code_folder)