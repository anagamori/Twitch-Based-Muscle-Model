%==========================================================================
% new_model_v11_parameterSearch.m
% Author: Akira Nagamori
% Last update: 2/22/119
% Descriptions
%   Inspired by Williams et al. 1998
%==========================================================================

close all
clear all
clc

%% Simulation parameters
Fs = 1000; %sampling frequency
T = 1/Fs;
time = 0:1/Fs:5; %simulation time
Lce = 1;

%% Parameters to be searched
target_CT = 90;


for m = 1501
    S = 7;
    C = 1;
    k_1 = 20;
    k_2 = 5;
    k_3 = 15;
    k_4 = 7;
    tau_1 = 0.005;
    tau_2 = 0.04;
    n = 1.8;
    k = 0.04;
    alpha = 4;
    
    r = randi([1 11],1,1);
    %     if r == 1
    %         S = 7+ (-S*0.8 + (S*0.8+S*0.8)*rand(1,1));
    %     elseif r == 2
    %         C = (1 + (4-1)*rand(1,1));
    %     elseif r == 3
    %         k_1 = 20 + (-k_1*0.8 + (k_1*0.8+k_1*0.8)*rand(1,1));
    %     elseif r == 4
    %         k_2 = 5 + (-k_2*0.8 + (k_2*0.8+k_2*0.8)*rand(1,1));
    %     elseif r == 5
    %         k_3 = 15 + (-k_3*0.8 + (k_3*0.8+k_3*0.8)*rand(1,1));
    %     elseif r == 6
    %         k_4 = 7 + (-k_4*0.8 + (k_4*0.8+k_4*0.8)*rand(1,1));
    %     elseif r == 7
    %         tau_1 = 0.005 + (-tau_1*0.8 + (tau_1*0.8+tau_1*0.8)*rand(1,1));
    %     elseif r == 8
    %         tau_2 = 0.04 + (-tau_2*0.8 + (tau_2*0.8+tau_2*0.8)*rand(1,1));
    %     elseif r == 9
    %         n = 1.8 + (-n*0.8 + (n*0.8+n*0.8)*rand(1,1));
    %     elseif r == 10
    %         k = 0.04 + (-k*0.8 + (k*0.8+k*0.8)*rand(1,1));
    %     elseif r == 11
    %         alpha = 4 + (-alpha*0.8 + (alpha*0.8+alpha*0.8)*rand(1,1));
    %     end
    
    S = 7+ (-S*0.8 + (S*0.8+S*0.8)*rand(1,1));
    
    C = (1 + (4-1)*rand(1,1));
    
    k_1 = 20 + (-k_1*0.8 + (k_1*0.8+k_1*0.8)*rand(1,1));
    
    k_2 = 5 + (-k_2*0.8 + (k_2*0.8+k_2*0.8)*rand(1,1));
    
    k_3 = 15 + (-k_3*0.8 + (k_3*0.8+k_3*0.8)*rand(1,1));
    
    k_4 = 7 + (-k_4*0.8 + (k_4*0.8+k_4*0.8)*rand(1,1));
    
    tau_1 = 0.005 + (-tau_1*0.8 + (tau_1*0.8+tau_1*0.8)*rand(1,1));
    
    tau_2 = 0.04 + (-tau_2*0.8 + (tau_2*0.8+tau_2*0.8)*rand(1,1));
    
    n = 1.8 + (-n*0.8 + (n*0.8+n*0.8)*rand(1,1));
    
    k = 0.04 + (-k*0.8 + (k*0.8+k*0.8)*rand(1,1));
    
    alpha = 4 + (-alpha*0.8 + (alpha*0.8+alpha*0.8)*rand(1,1));
    
    
    %% Run a twitch simulation and sweep simulation
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
            FR_test = [2 5 8 10 12 15 18 20 25 30 40 50 60 70 80 100 200]; %10:10:100];
        end
        
        %% initialization
        mean_exc = zeros(1,length(FR_test));
        p2p_exc = zeros(1,length(FR_test));
        
        %% Loop through test frequencies
        for f = 1:length(FR_test)
            FR = FR_test(f);
            %%  Generate spike train
            spike = zeros(1,length(time));
            
            if simulation_condition == 1
                % Generate a pulse
                spike(1*Fs) = 1;
            else
                % Generate spike train
                temp =spikeTrainGenerator(0:1/Fs:3,Fs,FR);
                spike(1*Fs:4*Fs) = temp;
            end
            
            
            %% Model parameters
            x = 0; % free calcium concentration
            y = 0; % concentraction of calcium bound to troponin
            z = 0; % porportion of bidning sites that formed cross-bridges (0-1)
            act = 0;
            
            %% Initialization
            x_vec = zeros(1,length(time));
            x_int_vec = zeros(1,length(time));
            y_vec = zeros(1,length(time));
            y_int_vec = zeros(1,length(time));
            z_vec = zeros(1,length(time));
            act_vec = zeros(1,length(time));
            
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
                    temp = conv(spike_temp,R_temp*(1+2*act^alpha));
                    %             R_i = 1 + (5-1)*exp(-(1/FR)/0.1);
                    %             temp = conv(spike_temp,R_temp*R_i);
                    R = R + temp(1:length(time));
                end
                
                x_dot = k_1*(C-x-y)*R(t) - k_2*x*((C*(S-1))+x+y)-(k_3*x+k_4*y)*(1-y);
                y_dot = (1-y)*(k_3*x-k_4*y);
                x = x_dot/Fs + x;
                y = y_dot/Fs + y;
                
                if y < 0
                    y_temp = 0;
                else
                    y_temp = y;
                end
                y_int = y_temp^n/(y_temp^n+k^n);
                
                %% Stage 3
                z_dot = (y_int-z)/tau_2; %(tau_temp/1000); %(1-z)*f_app - g_app*z;
                z = z_dot/Fs + z;
                
                %% Stage 4
                %
                act = z; % (z*(f_app_max+g_app)/f_app_max);
                
                %% Store variables
                %x_vec(t) = x(t);
                x_vec(t) = x;
                y_vec(t) = y;
                z_vec(t) = z;
                act_vec(t) = act;
                
            end
            
            mean_exc(f) = mean(act_vec(3.5*Fs:4*Fs));
            p2p_exc(f) = max(act_vec(3.5*Fs:4*Fs))-min(act_vec(3.5*Fs:4*Fs));
            
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
            %% Calculate twitch-tetanus ratio
            twitch2tetanus_ratio = p2p_exc(1)/mean_exc(f)
            %% Calculate the degree of fusion
            fusion = 1-p2p_exc/p2p_exc(1);
            
            %% Calculate FR_half
            FR_new = FR_test(1):0.1:FR_test(end);
            Af_new = spline(FR_test,mean_exc,FR_new);
            [~,loc] = min(abs(Af_new-0.5));
            FR_half = FR_new(loc);
            FR_half
            
            %% Calculate the desired activation-frequency relationship bassed on Song et al. (2008)
            f_eff = FR_new/FR_half;
            a_f = 0.56;
            n_f0 = 2.1;
            n_f1 = 5;
            n_f = n_f0 +n_f1* (1/Lce-1);
            Af_Song = 1-exp(-(f_eff./(a_f*n_f)).^n_f);
            
            %% Calculate error between the desired and generated activation-frequency relationship
            for j = 1:length(FR_test)
                [val,loc] = min(abs(FR_new-FR_test(j)));
                error_temp(j) = abs(Af_Song(loc)-mean_exc(j));
            end
            error = sum(error_temp)
            
            %% Plot results
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
            xlim([0 3])
            
            figure(8)
            plot(mean_exc./max(mean_exc),fusion,'LineWidth',1)
            xlabel('Activation','FontSize',14)
            ylabel('Fusion','FontSize',14)
            hold on
            plot(0:0.1:1,0:0.1:1,'--')
        end
        
    end
    
    Data = cell(2,12);
    Data(1,:) = {'t_0_100','t_100_50','t_40_10','P_t','t2t','FR_half','Error_CT','Error_Af','Test frequencies','Af','Fusion','Parameters'};
    Data{2,1} = t_0_100;
    Data{2,2} = t_100_50;
    Data{2,3} = t_40_10;
    Data{2,4} = pks;
    Data{2,5} = twitch2tetanus_ratio;
    Data{2,6} = FR_half;
    Data{2,7} = target_CT-t_0_100;
    Data{2,8} = error;
    Data{2,9} = FR_test/FR_half;
    Data{2,10} = mean_exc;
    Data{2,11} = fusion;
    Data{2,12} = [S,C,k_1,k_2,k_3,k_4,tau_1,tau_2,n,k,alpha];
    
    code_folder = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model';
    data_folder = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Data';
    
    cd(data_folder)
    save(['Data_' num2str(m)],'Data')
    cd(code_folder)
    
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
