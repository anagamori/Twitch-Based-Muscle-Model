close all
clear all
clc

dataFolder = '/Volumes/DATA2/Motor Unit Model Data/Twitch-based Muscle Model';
codeFolder = '/Users/akira/Documents/Github/Twitch-Based-Muscle-Model';
amp_vec = 0.1:0.1:1;

for i = 1:length(amp_vec)
    
    Data = cell(1,10);
    tic
    parfor j = 1:10
        output = twitchBasedMuscleModel_FCR(amp_vec(i));
        Data{j} = output;
    end
    toc
    cd (dataFolder)
    save(['Data' num2str(i)],'Data')
    cd (codeFolder)
    
    i
    
end