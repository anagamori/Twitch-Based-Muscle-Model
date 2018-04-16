close all
clear all
clc

Fs = 1000;
time = 0:1/Fs:2;

force_vec = zeros(1,length(time));
spikeTrain = zeros(1,length(time));
spikeTrain(1*Fs) = 1;
force = 0;
T = 0.03;
for t = 1:length(time)
    if spikeTrain(t) >= force
        T = 0.03;
    else
        T = 0.1;
    end
    force_dot = (spikeTrain(t)-force)/T;
    force = force_dot/Fs + force;
    
    force_vec(t) = force;
    
end

figure(1)
subplot(2,1,1)
plot(time,spikeTrain)
subplot(2,1,2)
plot(time,force_vec)