%==========================================================================
%PDRvsUth_nonlinear.m
% Author: Akira Nagamori
% Last update: 6/27/20
% Descriptions:
%   Plot the relationship between excitation and discharge rate of
%   individual motor units
%   Used to generate figure (onion_skinVsAHP)
%==========================================================================
close all
clear all
clc

cd('/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Development_Code/Data');
load('modelParameter_v3')
cd('/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Code for Figures')

%%
x = [11.025  21.277 50.097 72.34 93.617];
DR_Connelly = [12.269 15.231 22.212 28.135 41.779];
error_Connelly = [14.596 19.356 28.452 35.856 50.029];
error_Connelly = error_Connelly - DR_Connelly;

%%

nMU = 200;

amp_vec = [0.106 0.29 0.62 0.78 0.93];
mean_Force = zeros(1,length(amp_vec));
mean_FR = zeros(nMU,length(amp_vec));
CoV_FR = zeros(nMU,length(amp_vec));
for i = 1:length(amp_vec)
    cd('/Volumes/DATA2/PLOS_CB_Data/withTendon/Data_Connely')
    load(['Data_' num2str(i) '_' num2str(1)])
    cd('/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Code for Figures')
    if i <= 2
        Fs = 10000;
        duration = 10*Fs;
    elseif i > 2 && i <=3
        Fs = 15000;
        duration = 10*Fs;
    else
        Fs = 20000;
        duration = 10*Fs;
    end
    Force = output.ForceTendon(end-duration+1:end);
    mean_Force(i) = mean(Force)/444.5844;
    for n = 1:nMU
        spike_time = find(output.spike_train(n,end-duration+1:end));
        ISI = diff(spike_time)/(Fs/1000);
        mean_FR(n,i) = mean(1./ISI*1000);
        CoV_FR(n,i) = std(ISI)/mean(ISI)*100;
    end
    
end

%%
figure(1)
% boxplot(DR_MU_all,'positions',[10 25 50 75 100])
% hold on
% scatter(10*ones(size(DR_MU_all(:,1))).*(1+(rand(size(DR_MU_all(:,1)))-0.5)/10),DR_MU_all(:,1),'r','filled')

errorbar(mean_Force*100,nanmean(mean_FR),nanstd(mean_FR),'LineWidth',1,'Color',[37  65 178]/255)
hold on
errorbar(x,DR_Connelly,error_Connelly,'LineWidth',1,'Color','k')
xlabel('Synaptic Input (%Maximum)','FontSize',12)
ylabel('Dsicharge Rate (Hz)','FontSize',12)
legend('New Model','Connelly et al. 1999')
set(gca,'TickDir','out');
set(gca,'box','off')
% scatter(U_th_vec,DR_vec,'filled','k','LineWidth',0.5)
% hold on
% plot(U_th_vec,DR_Calc,'r','LineWidth',1)
% xlim([-1 11])
% xticks([0 0.05 0.1 0.15 0.2]*100)
% xlabel('Recruitment Threshold (%Maximum)','FontSize',8)
% ylabel('Discharge Rate (Hz)','FontSize',8)
% set(gca,'TickDir','out');
% set(gca,'box','off')
% ax = gca;
% ax.FontSize = 6;
% fig = gcf;
% fig.PaperUnits = 'inches';
% fig.PaperPosition = [0 0 3.34 3.34];




