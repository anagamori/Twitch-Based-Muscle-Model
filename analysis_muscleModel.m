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
data_folder = '/Volumes/DATA2/PLOS_CB_Data/noTendon/Model_onion_skin_N_100_Ur_80';
code_folder = '/Users/akira/Documents/Github/Twitch-Based-Muscle-Model';
figure_folder = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Figures';

%% 
amp_vec = [0.025 0.05 0.1:0.1:1];
Fs = 10000;
nTrial = 10;

Force_mat = zeros(nTrial,15*Fs+1);
mean_Force = zeros(nTrial,length(amp_vec));
std_Force = zeros(nTrial,length(amp_vec));
cov_Force = zeros(nTrial,length(amp_vec));
pxx = zeros(nTrial,201);
mean_pxx = zeros(length(amp_vec),201);
%% 
for j = 1:length(amp_vec)
    
    for i = 1:nTrial
        cd(data_folder)
        load(['Data_' num2str(j) '_' num2str(i)])
        cd(code_folder)
        Force =  output.Force;
        mean_Force(i,j) = mean(Force(5*Fs+1:end));
        std_Force(i,j) = std(Force(5*Fs+1:end));
        cov_Force(i,j) =  std_Force(i,j)/mean_Force(i,j)*100;
        [pxx(i,:),f] = pwelch(Force(5*Fs+1:end)-mean(Force(5*Fs+1:end)),[],[],0:0.5:100,Fs,'power');
        
        figure(11)
        plot(Force)
        hold on
        
        Force_mat(i,:) = Force;
    end
    mean_pxx(j+2,:) = mean(pxx);
    
    cd(data_folder)
    save(['Force_mat_' num2str(j)],'Force_mat')
    cd(code_folder)
end


%%
%close all
mean_mean_Force = mean(mean_Force);
figure(1)
plot([0 amp_vec],[0 mean(mean_Force)]./mean_mean_Force(end),'LineWidth',2)
xlabel('Activation (%)','FontSize',14)
ylabel('Force (%MVC)','FontSize',14)
set(gca,'TickDir','out');
set(gca,'box','off')

figure(2)
plot(amp_vec,mean(std_Force),'LineWidth',2)
xlabel('Mean Force (%)','FontSize',14)
ylabel('SD (N)','FontSize',14)
set(gca,'TickDir','out');
set(gca,'box','off')

figure(3)
% plot(amp_vec,mean(cov_Force),'LineWidth',2)
% hold on 
plot(amp_vec,mean(cov_Force),'LineWidth',2)
xlabel('Mean Force (%)','FontSize',14)
ylabel('CoV (%)','FontSize',14)
set(gca,'TickDir','out');
set(gca,'box','off')

figure(4)
plot(f,mean_pxx([1 3 5 8],:),'LineWidth',2)
xlabel('Frequency (Hz)','FontSize',14)
ylabel('Power (N^2)','FontSize',14)
set(gca,'TickDir','out');
set(gca,'box','off')
legend('10% MVC','30% MVC','50% MVC','80% MVC')
