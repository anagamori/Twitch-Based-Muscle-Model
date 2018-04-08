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

for i = 1 %2:length(amp_vec)
    cd (dataFolder)
    load(['Data_' num2str(i)],'Data')
    cd (codeFolder)
    
    Force = zeros(10,5*Fs);
   
    
    for j = 1:10
        output = Data{j};
        %Force(j,:) = sum(output.Force(:,5*Fs+2:end));   
        Force(j,:) = output.Tendon_Force(5*Fs+2:end);   
    end
    
    clear Data 
    
    cd (dataFolder)
    save(['Force_' num2str(i)],'Force')
    cd (codeFolder)
    i
    
    
end