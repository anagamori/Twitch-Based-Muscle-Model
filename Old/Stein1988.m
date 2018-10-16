close all
clear all
clc

%--------------------------------------------------------------------------
% a = 11.8;
% b_0 = 16.2;
% B = 12.8;
% k = 0.83;
% n = 3.7;
% b_1 = 0.66;
%
a = 9.4;
b_0 = 17.5;
B = 11.4;
k = 0.85;
n = 3.6;
b_1 = 0.6;


Fs = 1000;
time = 0:1/Fs:6;
t_temp = 0:1/Fs:3;


FR = [2]; % 5 10 15 20 30 40 50]; %round(3*FR_half(testingUnit)); %round(3*FR_half(testingUnit)) %round(2.7*FR_half(testingUnit)); % [2 5 10 15 20 25 30 35 40 45 48];

for f = 1:length(FR)
    FR_test = FR(f);
    %--------------------------------------------------------------------------
    % Generate spike train
    spikeTrain_temp = spikeTrainGenerator(t_temp,Fs,FR_test);
    spikeTrain = [zeros(1,1*Fs) spikeTrain_temp zeros(1,2*Fs)];
    %spikeTrain = zeros(1,length(time));
    %spikeTrain(1*Fs) = 1;
    spike_train_temp = zeros(1,length(spikeTrain));
    force = zeros(1,length(spikeTrain));
    force_2 = zeros(1,length(spikeTrain));
    x = zeros(1,length(spikeTrain));
    h = exp(-a*time);
    b_vec = zeros(1,length(spikeTrain));
    g_vec = zeros(1,length(spikeTrain));
    for t = 1:length(time)
        if t > 1
            b = b_0*(1-b_1*force(t-1)/B)^2;
        else
            b = b_0*(1-b_1*0/B)^2;
        end
        b_vec(t) = b;
        if spikeTrain(t) == 1
            spike_train_temp(t) = 1;
            q_temp = conv(spike_train_temp,h);
            q = q_temp(1:length(time));            
            x = q.^n./(q.^n+k^n);            
            g = exp(-b*time);
            force_temp = conv(x,g); %/(B*b)*exp(1);
            force = force_temp(1:length(time));
        end
        
        
    end
    
    
    figure(1)
    plot(time,force)
    hold on
    
    mean_Force(f) = mean(force(3*Fs:4*Fs));
    
end

function spikeTrain = spikeTrainGenerator(t,Fs,freq)

spikeTrain = zeros(1,length(t));
ISI = round(1/freq*Fs);
numSpikes = round(length(t)/ISI);
index = [1:numSpikes]*ISI;
index(index>length(t)) = [];
spikeTrain(index) = 1;
spikeTrain(1) = 1;

end
