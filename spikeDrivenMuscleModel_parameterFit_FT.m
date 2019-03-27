%==========================================================================
% spikeDriveMuscleModel_parameterFit_FT.m
% Author: Akira Nagamori
% Last update: 3/27/19
% Descriptions: 
%   Adjust model parameters to a unit with a target contraction time
%==========================================================================

close all
clear all
clc

%% Folder name
code_folder = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model';
data_folder = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Model Parameters/Model_1';

%% Simulation parameters
Fs = 5000; %sampling frequency
T = 1/Fs;
time = 0:1/Fs:5; %simulation time

Lce = 1; % Muscle length

%% Load a vector that contains contraction times of each unit
cd (data_folder)
load('CT')
cd (code_folder)

%%
first_MU = 221;
last_MU = 300;
Data_cell = cell(1,last_MU);

%% Weighting for optimization
weight_temp = ((70-5)*rand(1,last_MU)+5);
%parpool(10)

%% Test each unit
parfor j = first_MU:last_MU
    
    %% Contraction time of the unit to be optimized to
    target_CT = CT_sorted(j);
    
    %% Seed the initial parameter set
    param = [4,2,40,30,50,40,0.005,0.03,1.7,0.04,5];
    
    %% Run annealing curve algorithm for 6 iterations
    for k = 1:6
        rng shuffle
        Param_matrix = annealing_curve(param,k);
        
        r = randperm(11);
        
        %% Loop through all parameters
        for n = 1:11
            
            %% Loop through all perturbations
            error_long = zeros(1,3);
            for l = 1:3
                [S,C,k_1,k_2,k_3,k_4,tau_1,tau_2,N,K,alpha] = parameter_Assigning(param,Param_matrix,r(n),l);
                
                %% Run a twitch simulation and sweep simulation
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
                            
                            T = t_0_100/Fs;
                            P =  pks/T;
                            twitch_Milner_temp = P.*time.*exp(1-time/T);
                            twitch_Milner = conv(spike,twitch_Milner_temp);
                        end
                        
                    end
                    
                    if i == 2
                        %% Calculate twitch-tetanus ratio
                        twitch2tetanus_ratio = p2p_exc(1)/mean_exc(f);
                        %% Calculate the degree of fusion
                        fusion = 1-p2p_exc/p2p_exc(1);
                        
                        %% Calculate FR_half
                        FR_new = 0.1:0.1:FR_test(end);
                        Af_new = spline(FR_test,mean_exc,FR_new);
                        [~,loc] = min(abs(Af_new-0.5));
                        FR_half = FR_new(loc);
                        
                        %% Calculate the desired activation-frequency relationship bassed on Song et al. (2008)
                        f_eff = FR_new/FR_half;
                        a_f = 0.56;
                        n_f0 = 2.1;
                        n_f1 = 5;
                        n_f = n_f0 +n_f1* (1/Lce-1);
                        Af_Song = 1-exp(-(f_eff./(a_f*n_f)).^n_f);
                        
                        %% Calculate error between the desired and generated activation-frequency relationship
                        error_temp = error_calculation(Af_Song,Af_new,f_eff);
                        
                        weight = weight_temp(j);
                        error = sum(error_temp) + abs(target_CT-t_0_100) + weight*twitch2tetanus_ratio;
                        
                    end
                    
                end
                
                error_long(l) = error;
                [min_error,loc_min_error] = min(error_long);
            end
            [S,C,k_1,k_2,k_3,k_4,tau_1,tau_2,N,K,alpha] = parameter_Assigning(param,Param_matrix,r(n),loc_min_error);
            param =  [S,C,k_1,k_2,k_3,k_4,tau_1,tau_2,N,K,alpha];
        end
    end
    
    [Data_temp] = spikeDrivenMuscleModel_testFunction(param,Lce,0,'fast');
    Data_cell{j} = Data_temp
    
end

for MU_No = first_MU:last_MU
    Data = Data_cell{MU_No};
    cd(data_folder)
    save(['Data_' num2str(MU_No)],'Data')
    cd(code_folder)
    
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