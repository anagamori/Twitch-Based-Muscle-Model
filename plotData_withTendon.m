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
data_folder = '/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Data/withTendon/30_CoV';
code_folder = '/Users/akiranagamori/Documents/Github/Twitch-Based-Muscle-Model';
figure_folder = '/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Figures';

Fs = 20000;
amp_vec = 0.1:0.1:1;
time =0:1/Fs:15;
mean_Force = zeros(10,length(amp_vec));
std_Force = zeros(10,length(amp_vec));
cov_Force = zeros(10,length(amp_vec));
pxx = zeros(10,201);
mean_pxx = zeros(length(amp_vec),201);

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
    
    figure(11)
    plot(time,Force_mat)
    hold on
        
    mean_pxx(j,:) = mean(pxx);
    %clear Force_mat
end

%%
close all
mean_mean_Force = mean(mean_Force);
figure(1)
plot([0 amp_vec],[0 mean(mean_Force)]./mean_mean_Force(end),'LineWidth',2)
xlabel('Activation (%)','FontSize',14)
ylabel('Force (%MVC)','FontSize',14)
set(gca,'TickDir','out');
set(gca,'box','off')
cd (figure_folder)
saveas(gcf,'activation2meanForce_withTendon_30CoV','pdf')
cd (code_folder)

%%
figure(2)
hp1 = line(amp_vec,mean(std_Force),'LineWidth',2,'Color','b');
xlabel('Activation Level (%)','FontSize',14)
ylabel('SD (N)','FontSize',14)
ax1 = gca;
set(ax1,'XColor','b','YColor','b')
ax2 = axes('Position',get(ax1,'Position'), 'XAxisLocation','top','YAxisLocation','right', 'Color','none','XColor','r','YColor','r'); 
hp2 = line(mean(mean_Force)./mean_mean_Force(end),mean(std_Force),'LineWidth',2,'Color','r','Parent',ax2);
xlabel('Mean Force (%)','FontSize',14)
set(gca,'TickDir','out');
set(gca,'box','off')
cd (figure_folder)
saveas(gcf,'meanForce2SD_wtihTendon_30CoV','pdf')
cd (code_folder)

%%
figure(3)
hp3 = line(amp_vec,mean(cov_Force),'LineWidth',2,'Color','b');
xlabel('Activation Level (%)','FontSize',14)
ylabel('CoV(%)','FontSize',14)
ax3 = gca;
set(ax3,'XColor','b','YColor','b')
ax4 = axes('Position',get(ax3,'Position'), 'XAxisLocation','top','YAxisLocation','right', 'Color','none','XColor','r','YColor','r'); 
hp4 = line(mean(mean_Force)./mean_mean_Force(end),mean(cov_Force),'LineWidth',2,'Color','r','Parent',ax4);
xlabel('Mean Force (%)','FontSize',14)
set(gca,'TickDir','out');
set(gca,'box','off')
cd (figure_folder)
saveas(gcf,'meanForce2CoV_withTendon_30CoV','pdf')
cd (code_folder)

%%
figure(4)
plot(f,mean_pxx([1 3 5 8],:),'LineWidth',2)
xlabel('Frequency (Hz)','FontSize',14)
ylabel('Power (N^2)','FontSize',14)
set(gca,'TickDir','out');
set(gca,'box','off')
legend('10% MVC','30% MVC','50% MVC','80% MVC')
cd (figure_folder)
saveas(gcf,'pxx_withTendon_30CoV','pdf')
cd (code_folder)