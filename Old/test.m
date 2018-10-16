close all
clear all
clc

%--------------------------------------------------------------------------
a = 3000; %1/0.03; %9.4
b = 3; %1/0.05; %;17.5;
c = 400; %4*10^5;
d = 150;
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
FR = 10; %:2:50;

for f = 1 %:length(FR)
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
    A = 0;
    B = 0;
    C = 0;
    D = 0;
    force = 0;
    for t = 1:length(time)
        %C = spikeTrain(t) + exp(-a*T)*C;
        A = spikeTrain(t) + exp(-a*T)*A;
        A_vec(t) = A;
        B = A*T + exp(-b*T)*B;
        B_vec(t) = B;
        C = B*T + exp(-c*T)*C;
        C_vec(t) = C;
        D = C*T + exp(-d*T)*D;
        D_vec(t) = D;
    end
    
    
    figure(1)
    plot(time,D_vec)
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
