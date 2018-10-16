close all
clear all
clc

%--------------------------------------------------------------------------
A = 3.85*1000; %[N/ms]
R_0 = 2;
K_m = 0.19;
tau_c = 20/1000; %[ms]
tau_1 = 46.3/1000; %[ms]
tau_2 = 720/1000;
%--------------------------------------------------------------------------
Fs = 1000;
T = 1/Fs;
time = 0:1/Fs:6;
t_temp = 0:1/Fs:3;


%FR = [2 5 10 15 20 30 40 50]; %round(3*FR_half(testingUnit)); %round(3*FR_half(testingUnit)) %round(2.7*FR_half(testingUnit)); % [2 5 10 15 20 25 30 35 40 45 48];
FR = 10; %:1:50;

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
    C_N = zeros(1,length(time));
    F = 0;
    t_previous = 0;
    t_stamp = 0;
    for t = 1:length(time)
        if spikeTrain(t) == 1
            t_stamp = [t_stamp time(t)];
            for i = 2:length(t_stamp)
                if i == 2
                    R = 1;
                    C_N_temp = R*time/tau_c.*exp(-time/tau_c);
                else
                    R = 1 + (R_0-1)*exp(-(t_stamp(i)-t_stamp(i-1))/tau_c);
                    C_N_temp = R*time/tau_c.*exp(-time/tau_c);
                end
                C_N(t:end) = C_N(t:end) + C_N_temp(1:end-t+1);
            end
        end

        
        dF = A*C_N(t)/(K_m+C_N(t)) - F/(tau_1+tau_2*C_N(t)/(K_m+C_N(t)));
        F = dF*1/Fs + F;
        Force(t) = F;
    end
    
    
    figure(1)
    plot(time,C_N)
    
    figure(2)
    plot(time,Force)
    
    
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
