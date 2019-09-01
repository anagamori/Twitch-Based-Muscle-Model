%==========================================================================
% plotData_withTendon.m
% Author: Akira Nagamori
% Last update: 3/11/19
% Descriptions:
%==========================================================================
close all
clear all
clc

%%
condition = 'Model_8_20_CoV_50_Ur_Rec_3_PT_100';
data_folder = ['/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Data/withTendon/' condition];
code_folder = '/Users/akira/Documents/Github/Twitch-Based-Muscle-Model';
figure_folder = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Figures';

Fs = 10000;
amp_vec = [0.05 0.1:0.1:1];
time =0:1/Fs:15;
mean_Force = zeros(10,length(amp_vec));
std_Force = zeros(10,length(amp_vec));
cov_Force = zeros(10,length(amp_vec));
pxx = zeros(10,201);
mean_pxx = zeros(length(amp_vec),201);

trial_vec = 0:7; %0:10; %[0:6 8:10];
for k = 1:length(trial_vec)
    j = trial_vec(k);
    cd(data_folder)
    load(['Force_mat_' num2str(j)])
    cd(code_folder)
    for i = 1:10
        Force =  Force_mat(i,:);
        mean_Force(i,j+1) = mean(Force(5*Fs+1:end));
        std_Force(i,j+1) = std(Force(5*Fs+1:end));
        cov_Force(i,j+1) =  std_Force(i,j+1)/mean_Force(i,j+1)*100;
              
        [pxx(i,:),f] = pwelch(Force(5*Fs+1:end)-mean(Force(5*Fs+1:end)),rectwin(length(Force(5*Fs+1:end))),0,0:0.5:100,Fs,'power');
    end
    
    figure(11)
    plot(time,Force_mat)
    hold on
        
    mean_pxx(j+1,:) = mean(pxx);
    %clear Force_mat
end

cd(data_folder)
save('mean_Force','mean_Force')
save('std_Force','std_Force')
save('cov_Force','cov_Force')
save('mean_pxx','mean_pxx')
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
% cd (figure_folder)
% saveas(gcf,['activation2meanForce_withTendon_' condition],'pdf')
% cd (code_folder)

%%
figure(2)
errorbar(amp_vec,mean(std_Force),std(std_Force),'LineWidth',2,'Color','b');
xlabel('Activation Level/Mean Force (%)','FontSize',14)
ylabel('SD (N)','FontSize',14)
hold on 
errorbar(mean(mean_Force)./mean_mean_Force(end),mean(std_Force),std(std_Force),'LineWidth',2,'Color','r');
legend('Activation Level','Force Level','Location','northwest')
yticks([0.05 0.1 0.15 0.2 0.25])
set(gca,'TickDir','out');
set(gca,'box','off')
% cd (figure_folder)
% saveas(gcf,['meanForce2SD_wtihTendon_' condition],'pdf')
% cd (code_folder)

%%
figure(3)
errorbar(amp_vec,mean(cov_Force),std(cov_Force),'LineWidth',2,'Color','b');
xlabel('Activation Level/Mean Force (%)','FontSize',14)
ylabel('CoV (%)','FontSize',14)
hold on 
errorbar(mean(mean_Force)./mean_mean_Force(end),mean(cov_Force),std(cov_Force),'LineWidth',2,'Color','r');
legend('Activation Level','Force Level')
set(gca,'TickDir','out');
set(gca,'box','off')
% cd (figure_folder)
% saveas(gcf,['meanForce2CoV_withTendon_' condition],'pdf')
% cd (code_folder)

%%
figure(4)
plot(f,mean_pxx([1 3 5 8],:),'LineWidth',2)
xlabel('Frequency (Hz)','FontSize',14)
ylabel('Power (N^2)','FontSize',14)
set(gca,'TickDir','out');
set(gca,'box','off')
legend('10% MVC','30% MVC','50% MVC','80% MVC')
% cd (figure_folder)
% saveas(gcf,['pxx_withTendon_' condition],'pdf')
% cd (code_folder)