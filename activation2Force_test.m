%--------------------------------------------------------------------------
% activation2Force_test.m
% Author: Akira Nagamori
% Last update: 4/13/18
% Code descriptions
% Ojbective:
%   Find model parameters to realize activation to force relationship used
%   in Loeb's model (Song et al., 2008)
%   Firing rate vs. force
%   Firing rate vs. Fusion
%--------------------------------------------------------------------------

close all
clear all
clc
%--------------------------------------------------------------------------
% motor unit parameters
N_MU = 300; % number of motor units
i_MU = 1:N_MU; % index for motor units

Ur = 0.6; % recruitment threshold for the lastly recruited motor unit
Ur_1 = 0.01; % reruitment threshold for the first unit
f_RT = fit([1 N_MU]',[Ur_1 Ur]','exp1');
coeffs_f_RT = coeffvalues(f_RT);
U_th = coeffs_f_RT(1)*exp(coeffs_f_RT(2)*i_MU); % the resulting recruitment threshold for individual units

MFR1_MU = 8; %minimum firing rate of first unit
MFRn_MU = 14; %minimum firing rate of last unit
RTEn_MU = U_th(end)-Ur_1;  %recruitment threshold of last unit
MFR_MU = MFR1_MU + (MFRn_MU-MFR1_MU) * ((U_th-Ur_1)./RTEn_MU);
PFR_MU = 4*MFR_MU; %peak firing rate
FR_half = PFR_MU./2; % firing rate at which half of maximum tension is achieved
g_e = (PFR_MU(end) - MFR_MU(end))/(1-U_th(end));
index_fast = 239;

CT_n = 20; % contraction time of the largest motor unit
FR_half_n = FR_half(end); % f0.5 of the largest motor unit
CT = 1.5*(CT_n*FR_half_n)./FR_half; % contraction time of individual units
CT = CT - (CT(end)-CT_n);
CT = CT/1000;
RT = CT; % half-relaxation time

Fs = 1000;
time = 0:1/Fs:6;
t_temp = 0:1/Fs:3;
%--------------------------------------------------------------------------
% simulation parameters
Lce = 1;
Y = 1;
S = 0.96;
U_eff = 0;
alpha_eff = 0;

cv_MU = 0;
testingUnit =  1;
%FR = 2:4:60; %PFR_MU(testingUnit);
FR_test = [10 20 30 40 50 60]; %[2 5 10 15 20 25 30 35 40 45 50 55 60];

twitch2tetanus = 5;
force_half = twitch2tetanus/2;


%% Obtain non-corrected activation-force relationship
% twitch amplitude = 1
for k = 1:length(FR_test)
    n = testingUnit;
    
    amp = - MFR_MU(n)/g_e + U_th(n) + FR_test(k)/g_e;
    U = [zeros(1,1*Fs) amp*ones(1,3*Fs) zeros(1,2*Fs+1)];
    
    FR_mat = zeros(1,length(time));
    spike_train = zeros(1,length(time));
    force = zeros(1,length(time));
    for t = 1:length(time)
        if U(t) >= U_eff
            T = 0.03;
        else
            T = 0.15;
        end
        U_eff_dot = (U(t)-U_eff)/T;
        U_eff = U_eff_dot*1/Fs + U_eff;
        
        if t > 1
            FR = g_e.*(U_eff - U_th(n)) + MFR_MU(n);
            f_env = FR/FR_half(n);
            if FR < MFR_MU(n)
                FR = 0;
            elseif FR > PFR_MU(n)
                FR = PFR_MU(n);
            end
            
            T_alpha = 0.05*f_env;
            alpha = 1-(1-exp(-force(t)/force_half.^3));
            alpha_eff_dot = (alpha-alpha_eff)/T_alpha;
            alpha_eff = alpha_eff_dot/Fs + alpha_eff;
            noise_FR = FR;
            FR_mat(t) = FR;
            
%             if n <= index_fast
%                 Af = Af_slow_function(f_env,Lce,Y);
%                 Af_cor = Af_slow_correction_function(f_env,Lce,Y);
%             else
%                 Af = Af_fast_function(f_env,Lce,Y);
%                 Af_cor = Af_fast_correction_function(f_env,Lce,Y);
%             end
            
            spike_train_temp = zeros(1,length(time));
            if FR >= MFR_MU(n)
                if ~any(spike_train) % initial time
                    spike_train(t) = 1;
                    spike_train_temp(t) = 1;
                    mu = 1/FR;
                    Z = randn(1);
                    Z(Z>3.9) = 3.9;
                    Z(Z<-3.9) = -3.9;
                    noise = 1/noise_FR*cv_MU*Z;
                    spike_time_temp = (mu + noise)*Fs;
                    if spike_time_temp < 2*1000/Fs
                        spike_time_temp = 2;
                    end
                    spike_time = round(spike_time_temp) + t;
                    [twitch,~,~] = twitch_function(1,f_env,Lce,CT(n),RT(n),Fs);
                    force_temp = conv(spike_train_temp,twitch);
                    force = force + force_temp(1:length(time));
                    spike_time_previous = t;
                else
                    if spike_time == t
                        spike_train(t) = 1;
                        spike_train_temp(t) = 1;
                        mu = 1/FR;
                        Z = randn(1);
                        Z(Z>3.9) = 3.9;
                        Z(Z<-3.9) = -3.9;
                        noise = 1/noise_FR*cv_MU*Z;
                        spike_time_temp = (mu + noise)*Fs;
                        if spike_time_temp < 2*1000/Fs
                            spike_time_temp = 2;
                        end
                        spike_time = round(spike_time_temp) + t;
                        
                        ISI = (t-spike_time_previous)/Fs;
                        FR_temp = 1/ISI;
%                         if n <= index_fast
%                             Af = Af_slow_function(FR_temp/FR_half(n),Lce,Y);
%                             Af_cor = Af_slow_correction_function(FR_temp/FR_half(n),Lce,Y);
%                         else
%                             Af = Af_fast_function(FR_temp/FR_half(n),Lce,Y);
%                             Af_cor = Af_fast_correction_function(FR_temp/FR_half(n),Lce,Y);
%                         end
                        %ISI = mu;
                        [twitch,~,~] = twitch_function(alpha_eff,f_env,Lce,CT(n),RT(n),Fs);
                        force_temp = conv(spike_train_temp,twitch);
                        force = force + force_temp(1:length(time));
                        spike_time_previous = t;
                    elseif FR_mat(t-1) == 0
                        spike_train(t) = 1;
                        spike_train_temp(t) = 1;
                        
                        mu = 1/FR;
                        Z = randn(1);
                        Z(Z>3.9) = 3.9;
                        Z(Z<-3.9) = -3.9;
                        noise = 1/noise_FR*cv_MU*Z;
                        spike_time_temp = (mu + noise)*Fs;
                        if spike_time_temp < 2*1000/Fs
                            spike_time_temp = 2;
                        end
                        spike_time = round(spike_time_temp) + t;
                        [twitch,~,~] = twitch_function(1,f_env,Lce,CT(n),RT(n),Fs);
                        force_temp = conv(spike_train_temp,twitch);
                        force = force + force_temp(1:length(time));                       
                        spike_time_previous = t;
                    end
                end
            end
        end
        alpha_eff_vec(t) = alpha_eff;
    end
    
    figure(1)
    plot(time,force)
    hold on
    
    
    mean_Force(k) = mean(force(3*Fs+1:4*Fs));
    P2P_Force(k) = max(force(3*Fs+1:4*Fs)) - min(force(3*Fs+1:4*Fs));
    
    index_spike = find(spike_train(2*Fs+1:3*Fs)==1);
    index_spike_diff = diff(index_spike)./Fs;
    mean_ISI = mean(index_spike_diff);
    SD_ISI = std(index_spike_diff);
    CoV_ISI = SD_ISI/mean_ISI;
    mean_FR = mean(1./index_spike_diff);
    
    
end

P2P_Force = 1-P2P_Force/1;
% 
figure(3)
plot(FR_test,mean_Force/mean_Force(end)*100,'LineWidth',2)
% hold on
% plot([MFR_MU(testingUnit) MFR_MU(testingUnit)],[0 105],'b')
% plot([PFR_MU(testingUnit) PFR_MU(testingUnit)],[0 105],'b')
% ylim([0 105])
% % /FR_half(testingUnit)
% % hold on
% % plot(FR/FR_half(testingUnit),meanForce_2/meanForce_2(end))
xlabel('Frequency (Hz)','FontSize',14)
ylabel('Force (%)','FontSize',14)
% %legend('Force 1','Force 2')
% 
figure(4)
plot(FR_test,P2P_Force*100,'LineWidth',2)
% hold on
% plot([MFR_MU(testingUnit) MFR_MU(testingUnit)],[0 105],'b')
% plot([PFR_MU(testingUnit) PFR_MU(testingUnit)],[0 105],'b')
% ylim([0 105])
% % hold on
% % plot(FR,P2PForce_2,'LineWidth',2)
% xlabel('Frequency (Hz)','FontSize',14)
% ylabel('Degree of Fusion (%)','FontSize',14)
% % legend('Force 1','Force 2')
% 

%% function used
function [twitch,T1,T2_temp] = twitch_function(alpha,f_eff,Lce,CT,RT,Fs)
T1 = CT*Lce^2+CT*f_eff;
T2_temp = (RT + RT*f_eff)/Lce;
T2 = T2_temp/1.68;
t_twitch = 0:1/Fs:3;
f_1 = t_twitch./T1.*exp(1-t_twitch./T1);
f_2 = t_twitch./T2.*exp(1-t_twitch./T2);

twitch_temp = [f_1(1:round(T1*Fs+1)) f_2(round(T2*Fs+1):end)];
twitch = alpha*twitch_temp;

% [~,loc] = min(abs(twitch_temp(1:round(T1*Fs+1))-f0));
% twitch = twitch_temp(loc:end);
if length(twitch) < length(t_twitch)
    twitch = [twitch zeros(1,length(t_twitch)-length(twitch))];
else
    twitch = twitch(1:length(t_twitch));
end
end

function Af = Af_slow_function(f_eff,L,Y)
a_f = 0.56;
n_f0 = 2.1;
n_f1 = 5;
n_f = n_f0 + n_f1*(1/L-1);
Af = 1 - exp(-(Y*f_eff/(a_f*n_f))^n_f);
end



function Af = Af_fast_function(f_eff,L,S)
a_f = 0.56;
n_f0 = 2.1;
n_f1 = 3.3;
n_f = n_f0 + n_f1*(1/L-1);
Af = 1 - exp(-(S*f_eff/(a_f*n_f))^n_f);
end

function FF = frequency2Force_slow_function_new(f_env,L,Y)
a_f = 0.56;
n_f0 = 2.1; %2.33
n_f1 = 5;
n_f = n_f0 + n_f1*(1/L-1);
FF = 1 - exp(-(Y*f_env/(a_f*n_f))^n_f);
offset = 0.375*L - 0.1775;
alpha = 1-offset;
FF = FF*alpha+offset;
FF = FF/f_env;
if FF > 1
    FF = 1;
end
end

function FF = frequency2Force_fast_function_new(f_env,L,S)
a_f = 0.56;
n_f0 = 2.1;
n_f1 = 3.3;
n_f = n_f0 + n_f1*(1/L-1);
FF = 1 - exp(-(S*f_env/(a_f*n_f))^n_f);
offset = 0.375*L - 0.1775;
alpha = 1-offset;
FF = FF*alpha+offset;
FF = FF/f_env;
if FF > 1
    FF = 1;
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