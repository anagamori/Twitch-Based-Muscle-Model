%==========================================================================
% spikeDrivenMuscleModel_test.m
% Author: Akira Nagamori
% Last update: 3/27/19
% Descriptions
%   A new model driven by spike trains inspired by Williams et al. 1998
%   Generatetwitch response and test activation-frequency response with a
%   given parameter set
%==========================================================================

close all
clear all
clc

%%
code_folder = '/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model';
data_folder = '/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Model Parameters/Model_1';

%% 
test_folder = '/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Model Parameters/Model_1/ST';
MU_No = 1;
%% Simulation parameters
% Simulation parameters
Fs = 1000; %sampling frequency
h = 1/Fs;
time = 0:1/Fs:5; %simulation time

Lce = 1; % muscle length
%% Model parameters
cd(test_folder)
load(['MU_' num2str(MU_No)])
cd(code_folder)

S = parameter(1); %7;
C = parameter(2); %1.025;
k_1 = parameter(3); %14.625;
k_2 = parameter(4); %4.9375;
k_3 = parameter(5)*Lce + parameter(6); %17.41*Lce - 2.85;
k_4 = parameter(7)*Lce + parameter(8); %-7.67*Lce + 14.92;
tau_1 = parameter(9); %0.0051;
tau_2 = parameter(10); % 0.04;
N = parameter(11)*Lce + parameter(12); %-2.26*Lce + 4.20;
K = parameter(13)*Lce + parameter(14); %-0.044*Lce + 0.080;
alpha = parameter(15); %4.475;
%%
for i = 1:2
    if i == 1
        % Generate a pulse to record twitch response
        FR_test = 1;
    elseif i == 2
        % Generate a set of spike trains at multiple frequencies
        FR_test = [2 5 8 10 12 15 18 20 25 30 40 50 60 70 80 100 200]; %10:10:100];
    end
    %% initialization
    mean_exc = zeros(1,length(FR_test));
    p2p_exc = zeros(1,length(FR_test));
    
    %% Test each stimulus frequency
    for f = 1:length(FR_test)
        %% Generate spike train
        FR = FR_test(f);
        spike = zeros(1,length(time));
        if i == 1
            % Generate a pulse
            spike(1*Fs) = 1;
        else
            % Generate spike train
            temp =spikeTrainGenerator(0:1/Fs:3,Fs,FR);
            spike(1*Fs:4*Fs) = temp;
        end
        
        %%  initialization
        c = 0; % free calcium concentration
        cf = 0; % concentraction of calcium bound to troponin
        A = 0; % muscle activation
        
        R_vec = zeros(1,length(time));
        c_vec = zeros(1,length(time));
        x_int_vec = zeros(1,length(time));
        cf_vec = zeros(1,length(time));
        A_tilda_vec = zeros(1,length(time));
        A_vec = zeros(1,length(time));
        
        spike_temp = zeros(1,length(time));
        R_temp = exp(-time/tau_1);
        R = zeros(1,length(time));
        for t = 1:length(time)
            %% Stage 1
            % Calcium diffusion to sarcoplasm
            spike_temp = zeros(1,length(time));
            if spike(t) == 1
                spike_temp(t) = 1;
                temp = conv(spike_temp,R_temp*(1+2*A^alpha));
                R = R + temp(1:length(time));
            end
            %R = spike(t) + exp(-h/tau_1)*R; %*(1+3*A^alpha);
            
            %%
            c_dot = k_1*(C-c-cf)*R(t) - k_2*c*(S-C+c+cf)-(k_3*c-k_4*cf)*(1-cf);
            cf_dot = (1-cf)*(k_3*c-k_4*cf);
            c = c_dot/Fs + c;
            cf = cf_dot/Fs + cf;
            
            %% Stage 2
            % Cooperativity and saturation
            if cf < 0
                cf_temp = 0;
            else
                cf_temp = cf;
            end
            A_tilda = cf_temp^N/(cf_temp^N+K^N);
            
            %% Stage 3
            % First-order dynamics to muscle activation, A
            %tau_2 = tau_2_0*(1-0.8*A_tilda)^2;
            A_dot = (A_tilda-A)/tau_2;
            A = A_dot/Fs + A;
            
            %% Store variables
            %x_vec(t) = x(t);
            %R_vec(t) = R;
            c_vec(t) = c;
            cf_vec(t) = cf;
            A_tilda_vec(t) = A_tilda;
            A_vec(t) = A;
            
            
        end
        
        mean_exc(f) = mean(A_vec(3.5*Fs:4*Fs));
        p2p_exc(f) = max(A_vec(3.5*Fs:4*Fs))-min(A_vec(3.5*Fs:4*Fs));
        
        if i == 1
            %------------------------------------------------------------------
            % Twitch analysis
            %------------------------------------------------------------------
            % Find the peak amplitude of twitch and its location
            [pks,locs_0_100] = max(A_vec);
            % Display peak amplitude value
            pks
            % Time it takes from zero force to the peak
            t_0_100 = (locs_0_100-1*Fs)*1000/Fs % convert it in the unit of ms
            
            % Find a half peak amplitude
            peak_half = pks/2;
            % Find time from peak to half maximum
            [~,t_100_50] = min(abs(A_vec(locs_0_100:end)-peak_half));
            
            % Find the location of t_100_50
            locs_100_50 = locs_0_100+t_100_50;
            % Time it takes from peak to half maximum
            t_100_50 = t_100_50*1000/Fs
            
            % Find value of 40% of peak amplitude
            peak_40 = pks*0.4;
            % Find value of 10% of peak amplitude
            peak_10 = pks*0.1;
            
            % Find time when force reaches 40% maximum after peak
            [~,t_40] = min(abs(A_vec(locs_0_100:end)-peak_40));
            % Find time when force reaches 10% maximum after peak
            [~,t_10] = min(abs(A_vec(locs_0_100:end)-peak_10));
            
            % Time it takes from 40% peak to 10% peak
            t_40_10 = (t_10-t_40)*1000/Fs
            
            % Locatino of t_40_10
            locs_40_10 = locs_0_100+t_10;
            
            T = t_0_100/1000;
            P =  pks/T;
            twitch_Milner_temp = P.*time.*exp(1-time/T);
            twitch_Milner = conv(spike,twitch_Milner_temp);
        end
        
        figure(i)
        plot(time,A_vec,'LineWidth',1)
        hold on
        if i == 1
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
    
    if i == 2
        
        twitch2tetanus_ratio = p2p_exc(1)/mean_exc(f)
        fusion = 1-p2p_exc/p2p_exc(1);
        
        FR_new = 0.1:0.1:FR_test(end);
        Af_new = spline(FR_test,mean_exc,FR_new);
        [~,loc] = min(abs(Af_new-0.5));
        FR_half = FR_new(loc);
        f_eff = FR_new/FR_half;
        
        a_f = 0.56;
        n_f0 = 2.1;
        n_f1 = 5;
        n_f = n_f0 +n_f1* (1/Lce-1);
        Af_Song = 1-exp(-(f_eff./(a_f*n_f)).^n_f);
        
        FR_half
        
        %% Calculate error between the desired and generated activation-frequency relationship
        for j = 1:length(find(f_eff<3.5))
            error_temp(j) = abs(Af_Song(j)-Af_new(j));
        end
        error = sum(error_temp)
        
        figure(6)
        plot(FR_test/FR_half,mean_exc,'LineWidth',1)
        xlabel('Frequency (Hz)','FontSize',14)
        ylabel('Activation','FontSize',14)
        hold on
        plot(f_eff,Af_new)
        plot(f_eff,Af_Song,'color','k')
        xlim([0 3])
        legend('New','New Fit','Song')
        
        figure(7)
        plot(FR_test/FR_half,fusion,'LineWidth',1)
        xlabel('Frequency (Hz)','FontSize',14)
        ylabel('Fusion','FontSize',14)
        hold on
        xlim([0 3])
        
        figure(8)
        plot(mean_exc./max(mean_exc),fusion,'LineWidth',1)
        xlabel('Activation','FontSize',14)
        ylabel('Fusion','FontSize',14)
        hold on
        plot(0:0.1:1,0:0.1:1,'--','color','k')
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
