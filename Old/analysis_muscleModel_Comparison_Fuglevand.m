%==========================================================================
% plotData_noTendon.m
% Author: Akira Nagamori
% Last update: 3/11/19
% Descriptions:
%==========================================================================
close all
clear all
clc

%%
code_folder = '/Users/akiranagamori/Documents/Github/Twitch-Based-Muscle-Model';
figure_folder = '/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Figures';

Fs = 1000;
amp_vec = 0.1:0.1:1;
time =0:1/Fs:15;
mean_Force = zeros(10,length(amp_vec));
std_Force = zeros(10,length(amp_vec));
cov_Force = zeros(10,length(amp_vec));
pxx = zeros(10,201);
mean_pxx = zeros(length(amp_vec),201);
for k = 1:2
    if k == 1
        data_folder = '/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Data/noTendon/10_CoV_50_Ur_Rec_2';
    else
        data_folder = '/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Data/Fuglevand/Model_1';
    end
for j = 1:10
    cd(data_folder)
    load(['Force_mat_' num2str(j)])
    cd(code_folder)
    for i = 1:10
        Force =  Force_mat(i,:);
        mean_Force(i,j) = mean(Force(5*Fs+1:end));
        std_Force(i,j) = std(Force(5*Fs+1:end));
        cov_Force(i,j) =  std_Force(i,j)/mean_Force(i,j)*100;
              
        [pxx(i,:),f] = pwelch(Force(5*Fs+1:end)-mean(Force(5*Fs+1:end)),gausswin(5*Fs),0.9*5*Fs,0:0.5:100,Fs,'power');
    end

    mean_pxx(j,:) = mean(pxx);
    %clear Force_mat
end

%%
mean_mean_Force = mean(mean_Force);
figure(1)
plot([0 amp_vec],[0 mean(mean_Force)]./mean_mean_Force(end),'LineWidth',2)
xlabel('Activation (%)','FontSize',14)
ylabel('Force (%MVC)','FontSize',14)
hold on
set(gca,'TickDir','out');
set(gca,'box','off')
cd (figure_folder)
saveas(gcf,'activation2Force_noTendon_Recruitment_1vs2','pdf')
cd (code_folder)

%%
figure(2)
errorbar(mean(mean_Force)./mean_mean_Force(end),mean(std_Force/mean_mean_Force(end)*100),std(std_Force/mean_mean_Force(end)),'LineWidth',2);
xlabel('Mean Force (%)','FontSize',14)
ylabel('SD (%MVC)','FontSize',14)
hold on 
%yticks([0.05 0.1 0.15 0.2 0.25])
set(gca,'TickDir','out');
set(gca,'box','off')
cd (figure_folder)
saveas(gcf,'meanForce2SD_noTendon_Recruitment_1vs2','pdf')
cd (code_folder)

%%
figure(3)
errorbar(mean(mean_Force)./mean_mean_Force(end),mean(cov_Force),std(cov_Force),'LineWidth',2);
xlabel('Mean Force (%)','FontSize',14)
ylabel('CoV (%)','FontSize',14)
hold on 
set(gca,'TickDir','out');
set(gca,'box','off')
cd (figure_folder)
saveas(gcf,'meanForce2CoV_noTendon_Recruitment_1vs2','pdf')
cd (code_folder)

%%
figure(k+3)
plot(f,mean_pxx([1 3 5 8],:),'LineWidth',2)
xlabel('Frequency (Hz)','FontSize',14)
ylabel('Power (N^2)','FontSize',14)
set(gca,'TickDir','out');
set(gca,'box','off')
legend('10% MVC','30% MVC','50% MVC','80% MVC')
cd (figure_folder)
saveas(gcf,'pxx_noTendon_Recruitment_1vs2','pdf')
cd (code_folder)

end

figure(1)
legend('New Model','Fuglevand Model','Location','northwest')

figure(2)
legend('New Model','Fuglevand Model','Location','northwest')

figure(3)
legend('New Model','Fuglevand Model')
