close all
clear all
clc

dataFolder = '/Users/akiranagamori/Documents/GitHub/SDN Data';
codeFolder = '/Users/akiranagamori/Documents/Github/Twitch-Based-Muscle-Model';
amp_vec = 0.1:0.1:1;

for i = 2:length(amp_vec)
    
    Data = cell(1,10);
    tic
    %parpool(2)
    for j = 1:10
        output = twitchBasedMuscleModel_FDI_noTendon(amp_vec(i));
        Data{j} = output;
    end
    toc
    cd (dataFolder)
    save(['Data_FDI_noTendon_' num2str(i)],'Data')
    cd (codeFolder)
    
    i
    
end