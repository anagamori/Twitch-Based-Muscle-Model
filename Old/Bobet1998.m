%close all
clear all
clc

%--------------------------------------------------------------------------
a = 9.4; %9.4
b_0 = 17.5; %;17.5;
B = 11.4;
k = 0.85;
n = 3.6;
b_1 = 0.6;
%
% a = 40;
% b_0 = 40;
% B = 25;
% k = 0.8;
% n = 3.7;
% b_1 = 0.4;
% 

Fs = 10000;
T = 1/Fs;
time = 0:1/Fs:6;
t_temp = 0:1/Fs:3;

%FR = [2 5 10 15 20 30 40 50]; %round(3*FR_half(testingUnit)); %round(3*FR_half(testingUnit)) %round(2.7*FR_half(testingUnit)); % [2 5 10 15 20 25 30 35 40 45 48];
FR = 2; %:2:50;

for f = 1:length(FR)
    FR_test = FR(f);
    %--------------------------------------------------------------------------
    % Generate spike train
    spikeTrain_temp = spikeTrainGenerator(t_temp,Fs,FR_test);
    spikeTrain = [zeros(1,1*Fs) spikeTrain_temp zeros(1,2*Fs)];
    %spikeTrain = zeros(1,length(time));
    %spikeTrain(1*Fs) = 1;
    spike_train_temp = zeros(1,length(spikeTrain));
    Force = zeros(1,length(spikeTrain));
    force_2 = zeros(1,length(spikeTrain));
    x = zeros(1,length(spikeTrain));
    h = exp(-a*time);
    b_vec = zeros(1,length(spikeTrain));
    g_vec = zeros(1,length(spikeTrain));
    C = 0;
    force = 0;
    for t = 1:length(time)
        %C = spikeTrain(t) + exp(-a*T)*C;
        C = spikeTrain(t) + (1-exp(-a*T))*exp(-a*T)*C;
        C_vec(t) = C;
        %Q = 1-exp(-(C*0.9).^1.8);
        Q = C^n/(C^n+k^n);
        Q_vec(t) = Q;
        b = b_0*(1-b_1*force/B)^2;
        b_vec(t) = b;
        force = Q*T*b + exp(-b*T)*force;
        force_vec(t) = force;
        %force = Q + exp(-b*T)*force;
        
        Force(t) = B*force;
    end
    
    
    figure(1)
    plot(time,Force)
    hold on
    
    mean_Force(f) = mean(Force(3*Fs:4*Fs));
    p2p_Force(f) = max(Force(3*Fs:4*Fs))-min(Force(3*Fs:4*Fs)); 
    
    maxForce = max(Force);
    spike_time = find(spikeTrain == 1);
    end_stimulation = spike_time(end);
    endForce = max(Force(end_stimulation-0.5*Fs:end_stimulation+0.5*Fs));
    [~,loc_rise] = min(abs((Force(1*Fs+1:1.1*Fs)-maxForce)));
    t_rise_time(f) = (time(1*Fs+1+loc_rise)-time(1*Fs+1));
    [~,loc_fall] = min(abs((Force(end_stimulation:end_stimulation+0.3*Fs)-endForce/2)));
    t_fall_time(f) = (time(end_stimulation+loc_fall)-time(end_stimulation ));
    
    
    
end

figure(2)
plot(C_vec)
hold on 

Tw2Tet = p2p_Force(1)/B

B_half = B/2;
[~,loc] = min(abs((mean_Force-B_half)));

FR_half_reference = 14.4;
FR_reference = [3 5 8 10 15 20 30 40 50 60 80]/FR_half_reference;
Force_reference = [6.34 7.3 13.7 22.7 55.1 71.5 87.2, 93.7 96 98.9 100];

figure(2)
plot(FR/FR(loc),mean_Force./B*100,'LineWidth',1)
xlim([0 3])
hold on
plot(FR_reference,Force_reference)

figure(3)
plot(FR/FR(loc),(1-p2p_Force./max(p2p_Force))*100,'LineWidth',1)
xlim([0 3])
hold on

figure(4)
plot(mean_Force./B*100,(1-p2p_Force./max(p2p_Force))*100)
hold on 
plot(mean_Force./B*100,mean_Force./B*100)

function spikeTrain = spikeTrainGenerator(t,Fs,freq)

spikeTrain = zeros(1,length(t));
ISI = round(1/freq*Fs);
numSpikes = round(length(t)/ISI);
index = [1:numSpikes]*ISI;
index(index>length(t)) = [];
spikeTrain(index) = 1;
spikeTrain(1) = 1;

end
