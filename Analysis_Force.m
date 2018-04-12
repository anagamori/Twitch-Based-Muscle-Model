close all
clear all
clc

dataFolder = '/Volumes/DATA2/Motor Unit Model Data/Twitch-based Muscle Model';
codeFolder = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model';
amp_vec = 0.1:0.1:1;

Fs = 40000;

mean_Force_all = zeros(1,length(amp_vec));
SD_Force_all = zeros(1,length(amp_vec));
SD_Force_SD_all = zeros(1,length(amp_vec));
CoV_Force_all = zeros(1,length(amp_vec));
CoV_Force_SD_all = zeros(1,length(amp_vec));

for i = 1:length(amp_vec)
    cd (dataFolder)
    load(['Force_FDI_0_1_1_10_' num2str(i)],'Force')
    cd (codeFolder)
    
    mean_Force = zeros(1,10);
    SD_Force = zeros(1,10);
    CoV_Force = zeros(1,10);
    
    for j = 1:10
        Force_temp = Force(j,:);        
        mean_Force(j) = mean(Force_temp);
        SD_Force(j) = std(Force_temp);
        CoV_Force(j) = SD_Force(j)/mean_Force(j);
    end
    
    mean_Force_all(i) = mean(mean_Force);
    SD_Force_all(i) = mean(SD_Force);
    SD_Force_SD_all(i) = std(SD_Force);
    CoV_Force_all(i) = mean(CoV_Force);
    CoV_Force_SD_all(i) = std(CoV_Force);
    
end

figure(1)
plot([0 amp_vec],[0 mean_Force_all/mean_Force_all(end)*100])

figure(2)
errorbar(amp_vec*100, SD_Force_all,SD_Force_SD_all,'-','Color',[0.078,0,0.831],'LineWidth',2,'Marker','o','MarkerFaceColor',[0.078,0,0.831])
xlabel('Excitation Level (%)','FontSize',14,'fontweight','bold')
ylabel('SD (N)','FontSize',14,'fontweight','bold')
xlim([0 105])

figure(3)
errorbar(mean_Force_all/mean_Force_all(end)*100, SD_Force_all,SD_Force_SD_all,'-','Color',[0.078,0,0.831],'LineWidth',2,'Marker','o','MarkerFaceColor',[0.078,0,0.831])
xlabel('Mean Force (%)','FontSize',14,'fontweight','bold')
ylabel('SD (N)','FontSize',14,'fontweight','bold')
xlim([0 105])

figure(4)
errorbar(amp_vec*100, CoV_Force_all*100,CoV_Force_SD_all*100,'-','Color',[0.078,0,0.831],'LineWidth',2,'Marker','o','MarkerFaceColor',[0.078,0,0.831])
xlabel('Excitation Level (%)','FontSize',14,'fontweight','bold')
ylabel('CoV (%)','FontSize',14,'fontweight','bold')
xlim([0 105])
