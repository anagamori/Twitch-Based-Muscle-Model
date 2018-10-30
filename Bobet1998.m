%close all
clear all
clc

%--------------------------------------------------------------------------
a = 20; %9.4
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
time = 0:1/Fs:3;
t_temp = 0:1/Fs:3;

%FR = [2 5 10 15 20 30 40 50]; %round(3*FR_half(testingUnit)); %round(3*FR_half(testingUnit)) %round(2.7*FR_half(testingUnit)); % [2 5 10 15 20 25 30 35 40 45 48];
FR = 1; %2:100;

for f = 1:length(FR)
    FR_test = FR(f);
    %--------------------------------------------------------------------------
    % Generate spike train
    %spikeTrain_temp = spikeTrainGenerator(t_temp,Fs,FR_test);
    %spikeTrain = [zeros(1,1*Fs) spikeTrain_temp zeros(1,2*Fs)];
    spikeTrain = zeros(1,length(time));
    spikeTrain(1*Fs) = 1;
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
        C = spikeTrain(t) + exp(-a*T)*C;
        %C = spikeTrain(t) + (1-exp(-a*T))*exp(-a*T)*C;
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
    
    
    [pks,locs_peak] = max(Force);
    CT = locs_peak-1*Fs;
    CT = CT/Fs
    
    peak_half = pks/2;
    [~,HRT] = min(abs(Force(locs_peak:end)-peak_half));
    locs_hrt = locs_peak+HRT;
    HRT = HRT/Fs
    
    figure(1)
    plot(time,Force)
    hold on
     
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
