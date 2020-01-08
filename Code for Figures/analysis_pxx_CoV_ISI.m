close all
clear all
clc

data_directory = '/Volumes/DATA2/New_Model/Fuglevand/N_120_CV_20';
code_directory = '/Users/akira/Documents/Github/Twitch-Based-Muscle-Model/Code for Figures';
%load ('Input')
Fs = 1000;
t = 0:1/Fs:15;

T = 0.09;
t_twitch = 0:1/Fs:1;
twitch = t_twitch./T.*exp(1-t_twitch./T);
[pxx_twitch,f] = pwelch(twitch,[],[],0:0.1:30,Fs);
pxx_twitch_norm = pxx_twitch./sum(pxx_twitch);

for k = 0 %0:10 %length(amp_temp)
    trialN = k; 
    % predefine model parameters
    
  
    cd (data_directory)
    load(['Trial_' num2str(trialN)])
    cd (code_directory)
    for i = 1 %:10
        output = Data{i};
        Force = output.Force(1,end-5*Fs+1:end);
        spike_train = output.SpikeTrain(1,end-5*Fs+1:end);
        
        Force_conv_temp = conv(output.SpikeTrain(1,:),twitch);
        Force_conv_temp2 = Force_conv_temp(1:length(output.Force(1,:)));
        Force_conv = Force_conv_temp2(end-5*Fs+1:end);
        
        [pxx,f] = pwelch(Force-mean(Force),[],[],0:0.1:30,Fs);
        pxx = pxx./sum(pxx);
        
        [pxx_spike,~] = pwelch(spike_train-mean(spike_train),[],[],0:0.1:30,Fs);
        pxx_spike = pxx_spike./sum(pxx_spike);
        
        [pxx_conv,~] = pwelch(Force_conv-mean(Force_conv),[],[],0:0.1:30,Fs);
        pxx_conv = pxx_conv./sum(pxx_conv);
        
    end
    
    figure(1)
    plot(f,pxx)
    hold on
    plot(f,pxx_conv)
    %plot(f,pxx_spike.*pxx_twitch_norm./sum(pxx_spike.*pxx_twitch_norm))
    plot(f,pxx_spike)
%     hold on
%     
%     plot(f,pxx_conv)
end

