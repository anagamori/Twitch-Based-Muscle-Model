close all
clear all
clc

dataFolder = '/Users/akiranagamori/Documents/GitHub/SDN Data';
codeFolder = '/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model';
amp_vec = 0.1:0.1:1;

Fs = 1000;

mean_Force_all_1 = zeros(1,length(amp_vec));
SD_Force_all_1 = zeros(1,length(amp_vec));
SD_Force_SD_all_1 = zeros(1,length(amp_vec));
CoV_Force_all_1 = zeros(1,length(amp_vec));
CoV_Force_SD_all_1 = zeros(1,length(amp_vec));

mean_Force_all_2 = zeros(1,length(amp_vec));
SD_Force_all_2 = zeros(1,length(amp_vec));
SD_Force_SD_all_2 = zeros(1,length(amp_vec));
CoV_Force_all_2 = zeros(1,length(amp_vec));
CoV_Force_SD_all_2 = zeros(1,length(amp_vec));

for i = 1:length(amp_vec)
    cd (dataFolder)
    load(['Data_FDI_noTendon_' num2str(i)],'Data')  
    cd (codeFolder)
    
    Force = zeros(10,5*Fs);
    mean_Force_1 = zeros(1,10);
    SD_Force_1 = zeros(1,10);
    CoV_Force_1 = zeros(1,10);
       
    for j = 1:10
        output_temp = Data{j};  
        Force_temp = sum(output_temp.Force);
        Force(j,:) = Force_temp(5*Fs+2:end);
        mean_Force_1(j) = mean(Force(j,:));
        SD_Force_1(j) = std(Force(j,:));
        CoV_Force_1(j) = SD_Force_1(j)/mean_Force_1(j);        
        
    end
    
    mean_Force_all_1(i) = mean(mean_Force_1);
    SD_Force_all_1(i) = mean(SD_Force_1);
    SD_Force_SD_all_1(i) = std(SD_Force_1);
    CoV_Force_all_1(i) = mean(CoV_Force_1);
    CoV_Force_SD_all_1(i) = std(CoV_Force_1);
    
    cd (dataFolder)
    save(['Force_FDI_noTendon_' num2str(i)],'Force')
    cd (codeFolder)
    
end

%%
close 
figure(1)
plot([0 amp_vec]*100,[0 mean_Force_all_1]./mean_Force_all_1(end),'-','Color',[0.851,0.325,0.098],'LineWidth',2)
xlabel('Excitation Level (%)','FontSize',14,'fontweight','bold')
ylabel('Mean Force (%)','FontSize',14,'fontweight','bold')
%legend('w/ Tendon','w/o Tendon')

%%
close 
figure(2)
errorbar(amp_vec*100, SD_Force_all_1,SD_Force_SD_all_1,'-','Color',[0.851,0.325,0.098],'LineWidth',2,'Marker','o','MarkerFaceColor',[0.851,0.325,0.098])
xlabel('Excitation Level (%)','FontSize',14,'fontweight','bold')
ylabel('SD (N)','FontSize',14,'fontweight','bold')
xlim([0 105])
%legend('w/ Tendon','w/o Tendon')

%%
figure(3)
errorbar(mean_Force_all_1/mean_Force_all_1(end)*100, SD_Force_all_1,SD_Force_SD_all_1,'-','Color',[0.851,0.325,0.098],'LineWidth',2,'Marker','o','MarkerFaceColor',[0.851,0.325,0.098])
xlabel('Mean Force (%)','FontSize',14,'fontweight','bold')
ylabel('SD (N)','FontSize',14,'fontweight','bold')
xlim([0 105])
%legend('w/ Tendon','w/o Tendon')

%%
figure(4)
errorbar(mean_Force_all_1/mean_Force_all_1(end)*100, CoV_Force_all_1*100,CoV_Force_SD_all_1*100,'-','Color',[0.851,0.325,0.098],'LineWidth',2,'Marker','o','MarkerFaceColor',[0.851,0.325,0.098])
xlabel('Excitation Level (%)','FontSize',14,'fontweight','bold')
ylabel('CoV (%)','FontSize',14,'fontweight','bold')
xlim([0 105])
%legend('w/ Tendon','w/o Tendon')
