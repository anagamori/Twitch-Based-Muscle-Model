%==========================================================================
% analysis_muscleModel_individualUnit.m
% Author: Akira Nagamori
% Last update: 5/21/20
% Descriptions:
%==========================================================================
close all
clear all
clc

%%
condition = 'Model_11_var_CoV_80_Ur_Rec_3';
data_folder = ['/Volumes/DATA2/New_Model/withTendon/' condition];
data_folder_git = ['/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Data/withTendon/' condition];
code_folder = '/Users/akira/Documents/Github/Twitch-Based-Muscle-Model';
figure_folder = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Figures';

model_parameter_folder =  '/Users/akira/Documents/Github/Twitch-Based-Muscle-Model/Model Parameters/Model_11_Ur_05';

%%
cd(model_parameter_folder)
load('modelParameter')
cd(code_folder)

%%
amp_vec = [0.025 0.05 0.1:0.1:1];



nTrial = 1;
nMU = 200;
Force_mat = zeros(nMU,15*10000+1);
mean_Force = zeros(nMU,length(amp_vec));
std_Force = zeros(nMU,length(amp_vec));
cov_Force = zeros(nMU,length(amp_vec));
pxx = zeros(nMU,1001);
mean_pxx = zeros(length(amp_vec),1001);
mean_pxx_2 = zeros(length(amp_vec),1001);
mean_FR = zeros(200,length(amp_vec));
CoV_FR = zeros(200,length(amp_vec));
% %%
j_vec = [-1:8 10];
testUnit = 10;

for j = -1:10 %length(amp_vec)
    %j = j_vec(j+2);
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
    duration = 8*Fs;
    j
    tic
  

    for i = 1:nTrial
        cd(data_folder)
        load(['Data_' num2str(j) '_' num2str(i)])
        cd(code_folder)
        Force =  output.force;
        mean_Force(:,j+2) = mean(Force(:,end-duration+2:end),2);
        std_Force(:,j+2) = std(Force(:,end-duration+2:end),[],2);
        cov_Force(:,j+2) =  std_Force(:,j+2)./mean_Force(:,j+2)*100;
       
        for n = 1:200
            spike_time = find(output.spike_train(n,end-duration+1:end));
            ISI = diff(spike_time)/(Fs/1000);
            mean_FR(n,j+2) = mean(1./ISI*1000);
            CoV_FR(n,j+2) = std(ISI)/mean(ISI)*100;
        end
        
    end
    toc
 
end

%%

mean_mean_Force = mean(mean_Force);
force_plot = mean_Force(testUnit,:);
force_plot(force_plot<0.001) = NaN;
figure(1)
yyaxis left
plot([0 amp_vec]*100,[0 force_plot],'LineWidth',2)
ylabel('Force (N)','FontSize',14)
ylim([1 2.5])
yyaxis right 
plot([0 amp_vec]*100,[0 mean_FR(testUnit,:)],'LineWidth',2)
xlim([15 105])
ylim([8 18])
xlabel('Synaptic Input (%Maximum)','FontSize',14)
ylabel('Mean Discharge Rate (Hz)','FontSize',14)
set(gca,'TickDir','out');
set(gca,'box','off')

%%
std_plot = std_Force(testUnit,:);
std_plot(std_plot<0.001) = NaN;
figure(2)
plot(amp_vec*100,std_plot,'LineWidth',2)
%plot(amp_vec,std_Force,'LineWidth',2)
xlim([15 105])
xlabel('Synaptic Input (%Maximum)','FontSize',14)
ylabel('SD (N)','FontSize',14)
set(gca,'TickDir','out');
set(gca,'box','off')

%%
figure(3)
yyaxis left
plot(amp_vec*100,cov_Force(testUnit,:),'LineWidth',2)
xlim([15 105])
xlabel('Synaptic Input (%Maximum)','FontSize',14)
ylabel('CoV of Force(%)','FontSize',14)
yyaxis right 
plot(amp_vec*100,CoV_FR(testUnit,:),'LineWidth',2)
ylabel('CoV of ISIs(%)','FontSize',14)
set(gca,'TickDir','out');
set(gca,'box','off')

