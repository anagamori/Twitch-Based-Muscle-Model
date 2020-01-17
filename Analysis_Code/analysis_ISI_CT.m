close all
clear all
clc

%%
code_folder = '/Users/akira/Documents/Github/Twitch-Based-Muscle-Model/Analysis_Code';
model_parameter_folder =  '/Users/akira/Documents/Github/Twitch-Based-Muscle-Model/Model Parameters/Model_8';

%%
cd(model_parameter_folder)
load('modelParameter')
cd(code_folder)

FR_half = modelParameter.FR_half;
MDR = modelParameter.MDR;


Fs = 10000;
time = 0:1/Fs:15;

cv = 0; % coefficient of variation for interspike intervals

T = 0.07;
t_twitch = 0:1/Fs:1;
twitch = t_twitch./T.*exp(1-t_twitch./T);

spike_train = spike_train_generator(time,8,cv,Fs);

force = conv(spike_train,twitch);
force = force(5*Fs+1:length(spike_train));

mean_force = mean(force)
cov_force = std(force)/mean(force)

T_2 = 0.04;
t_twitch = 0:1/Fs:1;
twitch_2 = t_twitch./T_2.*exp(1-t_twitch./T_2);

spike_train = spike_train_generator(time,14,cv,Fs);

force_2 = conv(spike_train,twitch_2);
force_2 = force_2(5*Fs+1:length(spike_train));

mean_force_2 = mean(force_2)
cov_force_2 = std(force_2)/mean(force_2)


% spike_time = find(spike_train);
% ISI = diff(spike_time)/(Fs/1000);
% mean_FR = mean(1./ISI*1000)
% CoV_ISI = std(ISI)/mean(ISI)*100
%

figure(1)
plot(force)
hold on
plot(force_2)

function spike_train = spike_train_generator(time,DR,cv,Fs)

Z = randn(1,length(time));
Z(Z>3.9) = 3.9;
Z(Z<-3.9) = -3.9;

DR_vec = zeros(1,length(time));
spike_train = zeros(1,length(time));
spike_time = 0;

for t = 1:length(time)
    if t > 1
        DR_vec(t) = DR;
        index_1 = DR >= 0 & DR_vec(t-1) == 0;
        index_2 = DR >= 0 & spike_time ==t;
        index = or(index_1,index_2);
        
        if index
            
            if ~any(spike_train) % when the motor unit fires at the first time
                spike_train(t) = 1; % add a spike to the vector
                
                mu = 1/DR;
                spike_time_temp = (mu + mu*cv*Z(t))*Fs;
                if spike_time_temp <= 0.002*Fs
                    spike_time_temp = 0.002*Fs;
                end
                spike_time = round(spike_time_temp) + t;
                
            else % when the motor unit have already fired at least once
                if spike_time == t % when the motor unit fires
                    spike_train(t) = 1;
                    
                    
                    mu = 1/DR; % interspike interval
                    spike_time_temp = (mu + mu*cv*Z(t))*Fs; % interspike interval
                    if spike_time_temp <= 0.002*Fs
                        spike_time_temp = 0.002*Fs;
                    end
                    spike_time = round(spike_time_temp) + t;
                    
                elseif t > spike_time + round(1/DR*Fs)
                    spike_train(t) = 1;
                    
                    mu = 1/DR; % interspike interval
                    spike_time_temp = (mu + mu*cv*Z(t))*Fs; % interspike interval
                    if spike_time_temp <= 0.002*Fs
                        spike_time_temp = 0.002*Fs;
                    end
                    spike_time = round(spike_time_temp) + t;
                end
            end
        end
    end
end
end