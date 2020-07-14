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
data_folder = '/Volumes/DATA2/PLOS_CB_Data/noTendon/Model_singleUnit_constantCoV';
code_folder = '/Users/akira/Documents/Github/Twitch-Based-Muscle-Model';
figure_folder = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Figures';

%%
test_unit = 91;
%amp_vec = [0.025 0.05 0.1:0.1:1];
amp_vec = 0.01:0.01:1;
Fs = 10000;
nTrial = 10;

Force_mat = zeros(nTrial,5*Fs+1);
mean_Force = zeros(nTrial,length(amp_vec));
std_Force = zeros(nTrial,length(amp_vec));
cov_Force = zeros(nTrial,length(amp_vec));
std_Force_dt = zeros(nTrial,length(amp_vec));
cov_Force_dt = zeros(nTrial,length(amp_vec));
mean_FR = zeros(nTrial,length(amp_vec));
CoV_ISI = zeros(nTrial,length(amp_vec));
pxx = zeros(nTrial,201);
mean_pxx = zeros(length(amp_vec),201);
%%
for j = 1:length(amp_vec)
    
    for i = 1:nTrial
        cd(data_folder)
        load(['Data_' num2str(test_unit) '_' num2str(j) '_' num2str(i)])
        cd(code_folder)
        Force =  output.Force(5*Fs+1:end);
        
        window = 1*Fs;
        bp = [1:window:length(Force)];
        Force_dt = detrend(Force,1,bp);
        
        mean_Force(i,j) = mean(Force);
        std_Force(i,j) = std(Force);
        std_Force_dt(i,j) = std(Force_dt);
        cov_Force(i,j) =  std_Force(i,j)/mean_Force(i,j)*100;
        cov_Force_dt(i,j) =  std(Force_dt)/mean_Force(i,j)*100;
        [pxx(i,:),f] = pwelch(Force-mean(Force),[],[],0:0.5:100,Fs,'power');
        
        
        spike_time = find(output.spike_train(5*Fs+1:end));
        ISI = diff(spike_time)/(Fs/1000);
        mean_FR(i,j) = mean(1./ISI*1000);
        CoV_ISI(i,j) = std(ISI)/mean(ISI)*100;
        
        
        Force_mat(i,:) = Force;
    end
    mean_pxx(j,:) = mean(pxx);
    
    figure(11)
    plot(Force)
    hold on
    
    %     cd(data_folder)
    %     save(['Force_mat_' num2str(j)],'Force_mat')
    %     cd(code_folder)
end


cd(data_folder)
save(['mean_Force_' num2str(test_unit)],'mean_Force')
save(['std_Force_' num2str(test_unit)],'std_Force')
save(['std_Force_dt_' num2str(test_unit)],'std_Force_dt')
save(['cov_Force_' num2str(test_unit)],'cov_Force')
save(['cov_Force_dt_' num2str(test_unit)],'cov_Force_dt')
save(['mean_FR_' num2str(test_unit)],'mean_FR')
save(['CoV_ISI_' num2str(test_unit)],'CoV_ISI')
cd(code_folder)
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
