close all
clear all
clc

dataFolder = '/Volumes/DATA2/Motor Unit Model Data/Twitch-based Muscle Model';
codeFolder = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model';
amp_vec = 0.1:0.1:1;

%Fs = 1000;
Fs = 50000;

mean_Force_all = zeros(1,length(amp_vec));
SD_Force_all = zeros(1,length(amp_vec));
SD_Force_SD_all = zeros(1,length(amp_vec));
CoV_Force_all = zeros(1,length(amp_vec));
CoV_Force_SD_all = zeros(1,length(amp_vec));

for i = 2:length(amp_vec)
%     cd (dataFolder)
%     load(['Data_' num2str(i)],'Data')
%     cd (codeFolder)
    
    Force = zeros(10,5*Fs);
   
    cd (dataFolder)
    for j = 1:10
        load (['output_FCR_long_0_1_1_10_' num2str(i) '_' num2str(j) '.mat'])        
        Force(j,:) = sum(output.Force(:,5*Fs+2:end));   
        %Force(j,:) = output.Tendon_Force(5*Fs+2:end);   
        %Force_temp = output.Tendon_Force;
        Force_temp  = sum(output.Force(:,5*Fs+2:end));
        clear output
    end
    cd (codeFolder)
    
    figure(1)
    plot(Force_temp)
    hold on 
    
    cd (dataFolder)
    save(['Force_FCR_long_0_1_1_10_' num2str(i)],'Force')
    cd (codeFolder)
    i
    
    
end