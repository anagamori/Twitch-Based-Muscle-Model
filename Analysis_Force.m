close all
clear all
clc

dataFolder = '/Volumes/DATA2/Motor Unit Model Data/Twitch-based Muscle Model';
codeFolder = '/Users/akira/Documents/Github/Twitch-Based-Muscle-Model';
amp_vec = 0.1:0.1:1;

Fs = 1000;

mean_Force_all = zeros(1,length(amp_vec));
SD_Force_all = zeros(1,length(amp_vec));
SD_Force_SD_all = zeros(1,length(amp_vec));
CoV_Force_all = zeros(1,length(amp_vec));

for i = 1:length(amp_vec)
    cd (dataFolder)
    load(['Data' num2str(i)],'Data')
    cd (codeFolder)
    
    mean_Force = zeros(1,10);
    SD_Force = zeros(1,10);
    CoV_Force = zeros(1,10);
    
    for j = 1:10
        output = Data{j};
        Force = sum(output.Force(:,5*Fs+1:end));
        mean_Force(j) = mean(Force);
        SD_Force(j) = std(Force);
        CoV_Force(j) = SD_Force(j)/mean_Force(j);
    end
    
    mean_Force_all(i) = mean(mean_Force);
    SD_Force_all(i) = mean(SD_Force);
    SD_Force_SD_all(i) = std(SD_Force);
    CoV_Force_all(i) = mean(CoV_Force);
    
end

figure(1)
plot(amp_vec,mean_Force_all)

figure(2)
plot(amp_vec,SD_Force_all)