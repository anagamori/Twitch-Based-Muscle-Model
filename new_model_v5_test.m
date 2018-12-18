%==========================================================================
% new_model_v5_test.m
% Author: Akira Nagamori
% Last update: 12/18/18
% Descriptions: Two-stage model with a two-state cross-bridge model
%==========================================================================

close all
clear all
clc

%==========================================================================
% Simulation parameters
Fs = 1000; %sampling frequency
T = 1/Fs;
time = 0:1/Fs:4; %simulation time

simulation_condition = 1;

if simulation_condition == 1
    % Generate a pulse
    FR_test = 1;
elseif simulation_condition == 2
    % Generate a spike train at a given frequency
    FR_test = 20;
elseif simulation_condition == 3
    % Generate a set of spike trains at multiple frequencies
    FR_test = [2 5 10 15 20 30 50 100]; %10:10:100];
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
    x = zeros(1,length(time)); % proportion of Ca2+ bound to troponin (0-1)
    y = 0; % porportion of cross-bridge sites available for binding (0-1)
    z = 0; % porportion of bidning sites that formed cross-bridges (0-1)
    act = 0;
    
    f_app_max = 400;
    g_app = 10;
    
    tau_1 = 0.05; 
    tau_2 = 0.05;
    tau_3 = 0.03;
    tau_4 = 0.01;
    
    time_diffusion = time; % simulated time for changes in permeability
    diffusion1 = 1-exp(-time_diffusion./tau_1); %increase in permeability
    diffusion2 = exp(-time_diffusion./tau_2); %decrease in permeability
    diffusion1_conv = zeros(1,length(time));
    diffusion2_conv = zeros(1,length(time));

    n = 3.5;
    k = 0.2;
    
    a = 0.002755;
    b = 5.865;
    %==========================================================================
    % initialization
    spike_temp = zeros(1,length(time));
    spike_2_temp = zeros(1,length(time));
    x_int_vec = zeros(1,length(time));
    x_vec = zeros(1,length(time));
    y_vec = zeros(1,length(time));
    z_vec = zeros(1,length(time));
    act_vec = zeros(1,length(time));
    
    for t = 1:length(time)
        if spike(t) == 1
            spike_temp(t) = 1;            
            diffusion1_conv = conv(spike_temp, diffusion1);
            diffusion2_conv = conv(spike_temp, diffusion2);
            spike_temp = zeros(1,length(time));
            x_temp = diffusion1_conv.*diffusion2_conv;
        else
            x_temp = zeros(1,length(time));
        end
                  
        x = x + x_temp(1:length(time));
        x_int = x(t)^n/(x(t)^n+k^n);
%         beta = tau_3+0.1*act;
%         y_dot = (x_int-y)/beta; %beta; 
%         y = y_dot/Fs + y;
        y_int = a*exp(b*x_int);
        f_app = f_app_max*y_int+20*act;
%         z_dot = (y-z)/tau_4; 
%         z = z_dot/Fs + z;
        z_dot = (1-z)*f_app - g_app*z; 
        z = z_dot/Fs + z;
        act = z^2; 
        x_vec(t) = x(t);
        x_int_vec(t) = x_int;
        y_vec(t) = y;
        z_vec(t) = z;
        act_vec(t) = act;
        
    end
    
    mean_exc(f) = mean(act_vec(2*Fs:3*Fs));
    p2p_exc(f) = max(act_vec(2*Fs:3*Fs))-min(act_vec(2*Fs:3*Fs));
    
    if simulation_condition == 1
        %------------------------------------------------------------------
        % Twitch analysis        
        %------------------------------------------------------------------
        % Find the peak amplitude of twitch and its location
        [pks,locs_0_100] = max(act_vec);
        % Display peak amplitude value
        pks 
        % Time it takes from zero force to the peak 
        t_0_100 = (locs_0_100-1*Fs)*1000/Fs % convert it in the unit of ms
        
        % Find a half peak amplitude
        peak_half = pks/2;
        % Find time from peak to half maximum
        [~,t_100_50] = min(abs(act_vec(locs_0_100:end)-peak_half));
        
        % Find the location of t_100_50
        locs_100_50 = locs_0_100+t_100_50;
        % Time it takes from peak to half maximum 
        t_100_50 = t_100_50*1000/Fs
        
        % Find value of 40% of peak amplitude
        peak_40 = pks*0.4;
        % Find value of 10% of peak amplitude
        peak_10 = pks*0.1;
        
        % Find time when force reaches 40% maximum after peak
        [~,t_40] = min(abs(act_vec(locs_0_100:end)-peak_40));
        % Find time when force reaches 10% maximum after peak
        [~,t_10] = min(abs(act_vec(locs_0_100:end)-peak_10));
      
        % Time it takes from 40% peak to 10% peak
        t_40_10 = (t_10-t_40)*1000/Fs
        
        % Locatino of t_40_10
        locs_40_10 = locs_0_100+t_10;
        
    end
    
    figure(1)
    plot(time,act_vec,'LineWidth',1)
    hold on
    if simulation_condition == 1
        hold on
        plot([time(locs_0_100) time(locs_0_100)],[0 1],'--')
        text(time(locs_0_100),pks,num2str(t_0_100))
        hold on
        plot([time(locs_100_50) time(locs_100_50)],[0 1],'--')
        plot([time(locs_40_10) time(locs_40_10)],[0 1],'--')
        text(time(locs_100_50),peak_half,num2str(t_100_50))
        text(time(locs_40_10),peak_10,num2str(t_40_10))
        text(time(locs_0_100)-time(locs_0_100)*0.4,pks,num2str(pks))        
    end
    
end

figure(1)
xlabel('Time (s)','FontSize',14)
ylabel('Activation','FontSize',14)

if simulation_condition == 2
    figure(2)
    plot(time,x_vec,'LineWidth',1)
    xlabel('Time (s)','FontSize',14)
    ylabel('x','FontSize',14)
    
    figure(3)
    plot(time,x_int_vec,'LineWidth',1)
    xlabel('Time (s)','FontSize',14)
    ylabel('x_{int}','FontSize',14)
    
    figure(4)
    plot(time,y_vec,'LineWidth',1)
    xlabel('Time (s)','FontSize',14)
    ylabel('y','FontSize',14)

elseif simulation_condition == 3
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
    xlim([0 3])
    
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
