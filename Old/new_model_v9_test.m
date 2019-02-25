%==========================================================================
% new_model_v_test.m
% Author: Akira Nagamori
% Last update: 1/7/18
%   Modify the v8 to have a sequential bin
%
%==========================================================================

close all
clear all
clc

%==========================================================================
% Simulation parameters
Fs = 1000; %sampling frequency
T = 1/Fs;
time = 0:1/Fs:4; %simulation time

for i = 1:2
    if i == 1
        simulation_condition = 1;
    elseif i == 2
        simulation_condition = 3;
    end
    if simulation_condition == 1
        % Generate a pulse
        FR_test = 1;
    elseif simulation_condition == 2
        % Generate a spike train at a given frequency
        FR_test = 10;
    elseif simulation_condition == 3
        % Generate a set of spike trains at multiple frequencies
        FR_test = [2 5 8 10 12 15 18 20 25 30 50 100]; %10:10:100];
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
        % x = zeros(1,length(time)); % proportion of Ca2+ bound to troponin (0-1)
        x = 0;
        y = 0; % porportion of cross-bridge sites available for binding (0-1)
        z = 0; % porportion of bidning sites that formed cross-bridges (0-1)
        act = 0;
        
        % Stage 1
        tau_1 = 0.001;
        K_1 = 80;
        K_2_max = 150;
        
        n = 2.2;
        k = 0.05;
        
        alpha = 2;
        beta = 100;
        % Stage 3
        f_app_max = 20; % See eq. 3 and 4 and subsequent texts on Westerblad & Allen (1994)
        g_app = 40;
        
        %=========================================================================
        % initialization
        x_vec = zeros(1,length(time));
        x_int_vec = zeros(1,length(time));
        y_vec = zeros(1,length(time));
        y_int_vec = zeros(1,length(time));
        z_vec = zeros(1,length(time));
        act_vec = zeros(1,length(time));
        K_vec = zeros(1,length(time));
        spike_temp = zeros(1,length(time));
        R_temp = exp(-time/tau_1);
        R = zeros(1,length(time));
        
        for t = 1:length(time)
            %% Stage 1
            % Calcium diffusion to sarcoplasm
            spike_temp = zeros(1,length(time));
            if spike(t) == 1
                spike_temp(t) = 1;
                %temp = conv(spike_temp,R_temp);
                % R to depend on the normalized firing rate
                temp = conv(spike_temp,R_temp*(1+1*act^3));
                %             R_i = 1 + (5-1)*exp(-(1/FR)/0.1);
                %             temp = conv(spike_temp,R_temp*R_i);
                R = R + temp(1:length(time));
            end
            
            K_2 = K_2_max; %/(1+alpha*act);
            
            x_dot = K_1*R(t) - K_2*x;
            x = x_dot/Fs + x;
            
            x_int = x^n/(x^n+k^n);
            %         if x > 1
            %            x = 1;
            
            %% Stage 3
            % Dynamics of cross-bridge formation (Huxley, 1957)
%             f_app = f_app_max*x_int;
%             z_dot = (1-z)*f_app - g_app*z; 
            K = beta*x_int;
            z_dot  = -(g_app+f_app_max*(K/(K+1)))*z+f_app_max*(K/(K+1));
            z = z_dot/Fs + z;
            
            %% Stage 4
            %
            act = (z*(f_app_max+g_app)/f_app_max);           
            
            %% Store variables
            %x_vec(t) = x(t);
            x_vec(t) = x;
            x_int_vec(t) = x_int;
            y_vec(t) = y;
            z_vec(t) = z;
            act_vec(t) = act;
            K_vec (t) = K;
            
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
            
            T = t_0_100/Fs;
            P =  pks/T;
            twitch_Milner_temp = P.*time.*exp(1-time/T);
            twitch_Milner = conv(spike,twitch_Milner_temp);
        end
        
        figure(i)
        plot(time,act_vec,'LineWidth',1)
        hold on
        if simulation_condition == 1   
            plot(time,twitch_Milner(1:length(time)))
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
    
    figure(i)
    xlabel('Time (s)','FontSize',14)
    ylabel('Activation','FontSize',14)
    
    if simulation_condition == 2
        figure(3)
        plot(time,x_vec,'LineWidth',1)
        xlabel('Time (s)','FontSize',14)
        ylabel('x','FontSize',14)
        
        figure(4)
        plot(time,x_int_vec,'LineWidth',1)
        xlabel('Time (s)','FontSize',14)
        ylabel('x_{int}','FontSize',14)
        
        figure(5)
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
        FR_half
        
        figure(6)
        plot(FR_test/FR_half,mean_exc,'LineWidth',1)
        xlabel('Frequency (Hz)','FontSize',14)
        ylabel('Activation','FontSize',14)
        hold on
        plot(f_eff,Af_Song)
        xlim([0 3])
        
        figure(7)
        plot(FR_test/FR_half,fusion,'LineWidth',1)
        xlabel('Frequency (Hz)','FontSize',14)
        ylabel('Fusion','FontSize',14)
        hold on
        
        figure(8)
        plot(mean_exc./max(mean_exc),fusion,'LineWidth',1)
        xlabel('Activation','FontSize',14)
        ylabel('Fusion','FontSize',14)
        hold on
        plot(0:0.1:1,0:0.1:1,'--')
    end
    
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
