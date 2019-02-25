%==========================================================================
% new_model_v2_test_f2a.m
% Author: Akira Nagamori
% Last update: 10/31/18
%==========================================================================

%close all
clear all
clc

%==========================================================================
% Simulation parameters

Fs = 1000; %sampling frequency
T = 1/Fs;
time = 0:1/Fs:4; %simulation time

simulation_condition = 2;

if simulation_condition == 1
    % Generate a pulse
    FR_test = 1;
else
    FR_test = 10; %[2 5 10 15 20 30 50 100 200]; %10:10:100];
end


%==========================================================================
% initialization
mean_exc = zeros(1,length(FR_test));
p2p_exc = zeros(1,length(FR_test));

for f = 1:length(FR_test)
    FR = FR_test(f);
    %==========================================================================
    % Generate spike train
    spike = zeros(1,length(time));
    
    if simulation_condition == 1
        % Generate a pulse
        spike(1*Fs) = 1;
    else
        % Generate spike train
        temp =spikeTrainGenerator(0:1/Fs:2,Fs,FR);
        spike(1*Fs:3*Fs) = temp;
    end
    
    %==========================================================================
    % Model parameters
    x = 0; % proportion of Ca2+ bound to troponin (0-1)
    y = 0; % porportion of cross-bridge sites available for binding (0-1)
    z = 0; % porportion of bidning sites that formed cross-bridges (0-1)
    
    tau_act = 0.097;
    tau_deact = 0.133;
    
    %==========================================================================
    % initialization
    x_vec = zeros(1,length(time));
    
    
    for t = 1:length(time)
        if spike(t) >= x
           x_dot = (spike(t) - x)/tau_act;
        else
           x_dot = (spike(t) - x)/tau_deact; 
        end
        x = x_dot/Fs + x;
        x_vec(t) = x;
        
    end
    
    mean_exc(f) = mean(x_vec(2*Fs:3*Fs));
    p2p_exc(f) = max(x_vec(2*Fs:3*Fs))-min(x_vec(2*Fs:3*Fs));
    
    if simulation_condition == 1
        % Twitch analysis
        [pks,locs_peak] = max(x_vec);
        pks
        CT = (locs_peak-1*Fs)*1000/Fs
        
        peak_half = pks/2;
        [~,HRT] = min(abs(x_vec(locs_peak:end)-peak_half));
        locs_hrt = locs_peak+HRT;
        HRT = HRT*1000/Fs
    end
    
    figure(1)
    plot(time,x_vec)
    hold on
    if simulation_condition == 1
        hold on
        plot([time(locs_peak) time(locs_peak)],[0 1],'--')
        text(time(locs_peak),pks,num2str(CT))
        hold on
        plot([time(locs_hrt) time(locs_hrt)],[0 1],'--')
        text(time(locs_hrt),peak_half,num2str(HRT))
        text(time(locs_peak)-time(locs_peak)*0.4,pks,num2str(pks))
    end
    
end


if simulation_condition == 2
    twitch2tetanus_ratio = p2p_exc(1)/mean_exc(f)
    fusion = 1-p2p_exc/p2p_exc(1);
    
    FR_new = FR_test(1):0.1:FR_test(end);
    Af_new = spline(FR_test,mean_exc,FR_new);
    [~,loc] = min(abs(Af_new-0.5));
    FR_half = FR_new(loc);
    f_eff = FR_new/FR_half;
    Af_Song = 1-exp(-(f_eff./(0.56*2.1)).^2.1);
    
    figure(2)
    plot(FR_test/FR_half,mean_exc,'LineWidth',1)
    xlabel('Frequency (Hz)','FontSize',14)
    ylabel('Activation','FontSize',14)
    hold on
    plot(f_eff,Af_Song)
    %xlim([0 3])
    
    figure(3)
    plot(FR_test/FR_half,fusion,'LineWidth',1)
    xlabel('Frequency (Hz)','FontSize',14)
    ylabel('Fusion','FontSize',14)
    hold on
    
    figure(4)
    plot(mean_exc./max(mean_exc),fusion,'LineWidth',1)
    xlabel('Activation','FontSize',14)
    ylabel('Fusion','FontSize',14)
    hold on
    plot(0:0.1:1,0:0.1:1,'--')
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
