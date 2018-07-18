close all
clear all
clc

dataFolder = '/Volumes/DATA2/Motor Unit Model Data/Twitch-based Muscle Model';
codeFolder = '/Users/akira/Documents/Github/Twitch-Based-Muscle-Model';
amp_vec = 0.1:0.1:1;
Fs = 1000;
for i = 1:length(amp_vec)
    
    tic
    %parpool(2)
    for j = 1 %:10
        output = twitchBasedMuscleModel_FDI(amp_vec(i));
        
%         cd (dataFolder)
%         save(['output_FDI_noTendon_0_1_1_5_180_' num2str(i) '_' num2str(j)],'output')
%         cd (codeFolder)
%         
%         clear output
    end
    toc   
    
    figure()
    plot(output.time,sum(output.Force))
    hold on 
    plot([output.time(1) output.time(end)],[output.F0 output.F0])
   
    i
    
    output_force = sum(output.Force);
    meanForce(i) = mean(output_force(5*Fs:end));
    sdForce(i) = std(output_force(5*Fs:end));
    CoVForce(i)= sdForce(i)/meanForce(i);
    
end