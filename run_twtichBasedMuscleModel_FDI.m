close all
clear all
clc

dataFolder = '/Volumes/DATA2/Motor Unit Model Data/Twitch-based Muscle Model';
codeFolder = '/Users/akira/Documents/Github/Twitch-Based-Muscle-Model';
amp_vec = 0.1:0.1:1;

for i = 1:length(amp_vec)
    
    tic
    %parpool(2)
    for j = 1:10
        output = twitchBasedMuscleModel_FDI_noTendon(amp_vec(i));
        
        cd (dataFolder)
        save(['output_FDI_noTendon_0_1_1_5_180_' num2str(i) '_' num2str(j)],'output')
        cd (codeFolder)
        
        clear output
    end
    toc   
    
    i
    
end