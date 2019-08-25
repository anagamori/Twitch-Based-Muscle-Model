%==========================================================================
% test_noise_samplingRate.m
% Author: Akira Nagamori
% Last update: 8/25/19
% Descriptions:
%   Test how sampling rate affects variance, distribution, and power
%   spectrum of simulated noise.
%==========================================================================


close all
clear all
clc

for i = 1:2
    if i == 1
        Fs = 1000;
        time = 0:1/Fs:10;
        amp = 10000;
    else
        Fs = 10000;
        time = 0:1/Fs:10;
        amp = 100000;
    end
    
    x = 0;
    x_vec = zeros(1,length(time));
    
    
    for t = 1:length(time)
        [x] = noise(x,Fs,amp);
        x_vec(t) = x;
    end
   
    var(x_vec)
    
    [pxx,f] = pwelch(x_vec,[],[],0:0.5:500,Fs,'power');
    
    figure(1)
    histogram(x_vec,'Normalization','probability')
    hold on 
    
    figure(2)
    plot(f,pxx./sum(pxx))
    hold on
    
end

figure(1)
legend('1k Hz','10k Hz')

figure(2)
legend('1k Hz','10k Hz')


function [x] = noise(x,Fs,amp)
    D = amp;
    tau = 0.005; 
    chi = normrnd(0,1,[1,size(x,1)]);
    x_dot = -x/tau + sqrt(D)*chi;
    x = x_dot*1/Fs + x;
end