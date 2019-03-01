%==========================================================================
% new_model_parameterFit_ST_length.m
% Author: Akira Nagamori
% Last update: 2/22/119
%==========================================================================

close all
clear all
clc

%% Folder name
code_folder = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model';
data_folder = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Data/ST';
%% Simulation parameters
Fs = 1000; %sampling frequency
T = 1/Fs;
time = 0:1/Fs:5; %simulation time

%% Parameters to be searched
%Lce = 1.1;
for trialN = 41:100
    trialN
    cd(data_folder)
    load(['Data_' num2str(trialN)])
    cd(code_folder)
    
    FR_half_initial = Data{2,6};
    param_initial = Data{2,12};
    
    Data_temp = cell(1,40);
    parfor j = 1:40
        
        if j <= 10
            Lce = 0.8;
        elseif j > 10 && j <= 20
            Lce = 0.9;
        elseif j > 20 && j <= 30
            Lce = 1.1;
        elseif j > 30
            Lce = 1.2;
        end
        condition = 10+j;
        FR_half = FR_half_initial;
        param = param_initial; 
        for k = 1:6
            rng shuffle
            Param_matrix = annealing_curve(param,k);
            
            r = randperm(4);
            
            %% Loop through all parameters
            for n = 1:length(r)
                if r(n) == 1
                    index = 5;
                elseif r(n) == 2
                    index = 6;
                elseif r(n) == 3
                    index = 9;
                elseif r(n) == 4
                    index = 10;
                end
                %% Loop through all perturbations
                error_long = zeros(1,3);
                for l = 1:3
                    [S,C,k_1,k_2,k_3,k_4,tau_1,tau_2,N,K,alpha] = parameter_Assigning(param,Param_matrix,index,l);
                    
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
                                y_int = y_temp^N/(y_temp^N+K^N);
                                
                                %% Stage 3
                                z_dot = (y_int-z)/tau_2; %(tau_temp/1000); %(1-z)*f_app - g_app*z;
                                z = z_dot/Fs + z;
                                
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
                                
                                % Time it takes from zero force to the peak
                                t_0_100 = (locs_0_100-1*Fs)*1000/Fs; % convert it in the unit of ms
                                
                                % Find a half peak amplitude
                                peak_half = pks/2;
                                % Find time from peak to half maximum
                                [~,t_100_50] = min(abs(act_vec(locs_0_100:end)-peak_half));
                                
                                % Find the location of t_100_50
                                locs_100_50 = locs_0_100+t_100_50;
                                % Time it takes from peak to half maximum
                                t_100_50 = t_100_50*1000/Fs;
                                
                                % Find value of 40% of peak amplitude
                                peak_40 = pks*0.4;
                                % Find value of 10% of peak amplitude
                                peak_10 = pks*0.1;
                                
                                % Find time when force reaches 40% maximum after peak
                                [~,t_40] = min(abs(act_vec(locs_0_100:end)-peak_40));
                                % Find time when force reaches 10% maximum after peak
                                [~,t_10] = min(abs(act_vec(locs_0_100:end)-peak_10));
                                
                                % Time it takes from 40% peak to 10% peak
                                t_40_10 = (t_10-t_40)*1000/Fs;
                                
                                % Locatino of t_40_10
                                locs_40_10 = locs_0_100+t_10;
                                
                                T = t_0_100/Fs;
                                P =  pks/T;
                                twitch_Milner_temp = P.*time.*exp(1-time/T);
                                twitch_Milner = conv(spike,twitch_Milner_temp);
                                
                            elseif simulation_condition == 3
                                %% Calculate twitch-tetanus ratio
                                twitch2tetanus_ratio = p2p_exc(1)/mean_exc(f);
                                %% Calculate the degree of fusion
                                fusion = 1-p2p_exc/p2p_exc(1);
                                
                                %% Calculate FR_half
                                FR_new = 0.1:0.1:FR_test(end);
                                Af_new = spline(FR_test,mean_exc,FR_new);
                                
                                %% Calculate the desired activation-frequency relationship bassed on Song et al. (2008)
                                f_eff = FR_new/FR_half;
                                a_f = 0.56;
                                n_f0 = 2.1;
                                n_f1 = 5;
                                n_f = n_f0 +n_f1* (1/Lce-1);
                                Af_Song = 1-exp(-(f_eff./(a_f*n_f)).^n_f);
                                
                                %% Calculate error between the desired and generated activation-frequency relationship
                                error_temp = error_calculation(Af_Song,Af_new,f_eff);
                                error = sum(error_temp);
                            end
                            
                        end
                    end
                    error_long(l) = error;
                    [min_error,loc_min_error] = min(error_long);
                end
                [S,C,k_1,k_2,k_3,k_4,tau_1,tau_2,N,K,alpha] = parameter_Assigning(param,Param_matrix,index,loc_min_error);
                param =  [S,C,k_1,k_2,k_3,k_4,tau_1,tau_2,N,K,alpha];
            end
        end
        
        [Data_temp{j}] = model_test(param,Lce,FR_half,'slow');
        
    end
    for trial = 1:40
        Data  = Data_temp{trial};
        cd(data_folder)
        save(['Data_' num2str(trialN) '_' num2str(trial+10)],'Data')
        cd(code_folder)
    end
end
%%
function spikeTrain = spikeTrainGenerator(t,Fs,freq)

spikeTrain = zeros(1,length(t));
ISI = round(1/freq*Fs);
numSpikes = round(length(t)/ISI);
index = [1:numSpikes]*ISI;
index(index>length(t)) = [];
spikeTrain(index) = 1;
spikeTrain(1) = 1;

end

function Y = annealing_curve(x,k)

perturbation_amp = 0.2/2;
S = 7;
C = 1;
k_1 = 20;
k_2 = 5;
k_3 = 15;
k_4 = 7;
tau_1 = 0.005;
tau_2 = 0.04;
N = 1.8;
K = 0.04;
alpha = 4;

Y(1,1) = x(1) +  (S*perturbation_amp)./2.^(k-2);
Y(2,1) = x(1) -  (S*perturbation_amp)./2.^(k-2);
Y(3,1) = x(1);

Y(1,2) = x(2) +  (C*perturbation_amp)./2.^(k-2);
Y(2,2) = x(2) -  (C*perturbation_amp)./2.^(k-2);
Y(3,2) = x(2);

Y(1,3) = x(3) +  (k_1*perturbation_amp)./2.^(k-2);
Y(2,3) = x(3) -  (k_1*perturbation_amp)./2.^(k-2);
Y(3,3) = x(3);

Y(1,4) = x(4) +  (k_2*perturbation_amp)./2.^(k-2);
Y(2,4) = x(4) -  (k_2*perturbation_amp)./2.^(k-2);
Y(3,4) = x(4);

Y(1,5) = x(5) +  (k_3*perturbation_amp)./2.^(k-2);
Y(2,5) = x(5) -  (k_3*perturbation_amp)./2.^(k-2);
Y(3,5) = x(5);

Y(1,6) = x(6) +  (k_4*perturbation_amp)./2.^(k-2);
Y(2,6) = x(6) -  (k_4*perturbation_amp)./2.^(k-2);
Y(3,6) = x(6);

Y(1,7) = x(7) +  (tau_1*perturbation_amp)./2.^(k-2);
Y(2,7) = x(7) -  (tau_1*perturbation_amp)./2.^(k-2);
Y(3,7) = x(7);

Y(1,8) = x(8) +  (tau_2*perturbation_amp)./2.^(k-2);
Y(2,8) = x(8) -  (tau_2*perturbation_amp)./2.^(k-2);
Y(3,8) = x(8);

Y(1,9) = x(9) +  (N*perturbation_amp)./2.^(k-2);
Y(2,9) = x(9) -  (N*0.2)./2.^(k-2);
Y(3,9) = x(9);

Y(1,10) = x(10) +  (K*perturbation_amp)./2.^(k-2);
Y(2,10) = x(10) -  (K*perturbation_amp)./2.^(k-2);
Y(3,10) = x(10);

Y(1,11) = x(11) +  (alpha*perturbation_amp)./2.^(k-2);
Y(2,11) = x(11) -  (alpha*perturbation_amp)./2.^(k-2);
Y(3,11) = x(11);

end

function [var1,var2,var3,var4,var5,var6,var7,var8,var9,var10,var11] = parameter_Assigning(x,M,index_1,index_2)

x(index_1) = M(index_2,index_1);
var1 = x(1);
var2 = x(2);
var3 = x(3);
var4 = x(4);
var5 = x(5);
var6 = x(6);
var7 = x(7);
var8 = x(8);
var9 = x(9);
var10 = x(10);
var11 = x(11);
end

function error = error_calculation(vec_1,vec_2,reference)
for i = 1:length(find(reference<2.5))
    error(i) = abs(vec_1(i)-vec_2(i));
end
end