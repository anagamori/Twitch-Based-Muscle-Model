%--------------------------------------------------------------------------
% activation2force.m
% Author: Akira Nagamori
% Last update: 10/15/18
%--------------------------------------------------------------------------

close all
clear all
clc
%--------------------------------------------------------------------------
% 1) Determine the number of motor units in a pool
N_MU = 300; % number of motor units
i_MU = 1:N_MU; % index for motor units

%--------------------------------------------------------------------------
% 2) Determine the maximum force of muscle
L0 = 5.1;
density = 1.06; % muscle density [g/cm^3]
mass = 0.0055; % muscle mass [kg]
PCSA = (mass*1000)/(density*L0); % PCSA of muscle
sigma = 31.8; % specific tension
F0 = PCSA*sigma;

%--------------------------------------------------------------------------
% 3) Determine the distribution of peak tetanic tension of each motor unit
RP_MU = 25; %range of peak tetanic tension between the smallest and largest motor units [ratio to the smallest]
b_MU = log(RP_MU)/N_MU; %coefficient to establish the exponential distribution of peak tetanic tension
PTf_MU = exp(b_MU*i_MU); %peak tetanic tension of individual motor units [ratio to the smallest]
PTi = PTf_MU./sum(PTf_MU)*F0; %correct PTf_MU to the unit of Newton [Newton]

%--------------------------------------------------------------------------
% 4) Determine the index for the largest slow-twitch fibers in a pool
F_pcsa_slow = 0.5; % fractional PSCA of slow-twitch motor units [0-1]
[~, index_slow] = min(abs(cumsum(PTi) - F0*F_pcsa_slow)); % index for the largest motor unit slow-twitch fiber in a pool

%--------------------------------------------------------------------------
% 5) Determine twithc-tetanus ratio
twitch2tetanus_ratio = 0.3;

%--------------------------------------------------------------------------
% 6) Determine the recruitment threshold of motor units [0 1]
Ur = 0.6; % recruitment threshold for the last-recruited motor unit (i.e., the largest motor unit)
Ur_1 = 0.01; % reruitment threshold for the first unit
b_Ur = log(Ur/Ur_1)/N_MU; %coefficient to establish the exponential distribution of recruitment threshold
U_th = exp(b_Ur*i_MU)/100; % recruitment threshold for each motor unit [0 1]

%--------------------------------------------------------------------------
% 7) Determine the minimum and maximum firing of motor units
MFR1_MU = 8; %minimum firing rate of first unit
MFRn_MU = 14; %minimum firing rate of last unit
RTEn_MU = U_th(end)-Ur_1;  %difference in the recruitment threshold of the smallest and largest motor units
MFR_MU = MFR1_MU + (MFRn_MU-MFR1_MU) * ((U_th-Ur_1)./RTEn_MU); %minimum firing rate of all motor units as a linear function of recruitment thresholds
PFR_MU = 4*MFR_MU; %peak firing rate
FR_half = PFR_MU./2; % firing rate at which half of maximum tetanic tension is achieved

%--------------------------------------------------------------------------
% 8) Determine the contraction time and half-reflaxation time
% contraction time was derived from the linear relationship with the
% FR_half (Thomas et al., 1991, Kernell et al. 1983b)
CT_n = 30; %contraction time of first unit [ms]
CT_1 = 90; %contraction time of lasdt unit [ms]
slope_CT = (CT_1 - CT_n)/(FR_half(1)-FR_half(end)); %slope for the linear relationship between contraction time and FR_half
intercept_CT = CT_1-slope_CT*FR_half(1); %intercept for the linear relationship between contraction time and FR_half
CT = slope_CT*FR_half+intercept_CT; %contraction time of all motor units
CT = CT/1000; %convert the contraction time into second [s]
RT = CT; %half-relaxation time [s]

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Simulation
Fs = 1000; %sampling frequency
time = 0:1/Fs:5; %simulation time
time_active = 0:1/Fs:3; %time window in which a motor unit discharges
testingUnit =  1; %index for the testing motor unit
Lce = 1; % muscle length
Y = 1;
S = 1;

FR = [2 5:3:3*FR_half(testingUnit) 3*FR_half(testingUnit)];
Af_old = length(FR);

for i = 1:length(FR)
    %force = zeros(1,length(time));
    f_env = FR(i)/FR_half(testingUnit);
    spikeTrain_temp = spikeTrainGenerator(time_active,Fs,FR(i));
    spikeTrain = [zeros(1,1*Fs) spikeTrain_temp zeros(1,1*Fs)];
    
    f_env_2 = [zeros(1,1*Fs) f_env*ones(1,length(spikeTrain_temp)) zeros(1,1*Fs)];
    f_int = 0;
    f_eff = 0;
    f_eff_dot = 0;
    f_eff_vec = zeros(1,length(f_env_2));
    Af_vec = zeros(1,length(f_env_2));
    
    if testingUnit <= index_slow
        Af_old(i) = Af_slow_function(f_env,Lce,Y);
        Af = Af_slow_function(f_env,Lce,Y);
    else
        Af_old(i) = Af_fast_function(f_env,Lce,S);
        Af = Af_fast_function(f_env,Lce,S);
    end
    
    for t = 1:length(f_env_2)
        if f_eff_dot >= 0
            T_f = 0.0343 + 0.0227*f_env;
        else
            T_f = (0.047+0.0252*Af_old(i));
        end
        f_int_dot = (f_env_2(t) - f_int)/T_f;
        f_int = f_int_dot*1/Fs + f_int;
        f_eff_dot = (f_int - f_eff)/T_f;
        f_eff = f_eff_dot*1/Fs + f_eff;
        f_eff_vec(t) = f_eff;
        if testingUnit <= index_slow
            Af_temp = Af_slow_function(f_env,Lce,Y);
        else
            Af_temp = Af_fast_function(f_env,Lce,S);
        end
        Af_vec(t) = Af_temp;
    end
    
    a = 12.56*exp(-4.229*Lce); %2*rand(1); 1 for Lce = 0.6, 0.4 for Lce = 0.8, 0.2 for Lce = 1.2
    b = 12.56*exp(-4.229*Lce);
    
    T1 = CT(testingUnit); %+CT(testingUnit)*Af_old(i);% *1.1;
    T2_temp = RT(testingUnit);% + RT(testingUnit)*Af_old(i));
    
    T2 = T2_temp/1.68;
    t_twitch = 0:1/Fs:5;
    f_1 = t_twitch./T1.*exp(1-t_twitch./T1);
    f_2 = t_twitch./T2.*exp(1-t_twitch./T2);
    
    twitch = [f_1(1:round(T1*Fs+1)) f_2(round(T2*Fs+1):end)];
    if length(twitch) < length(t_twitch)
        twitch = [twitch zeros(1,length(t_twitch)-length(twitch))];
    else
        twitch = twitch(1:length(t_twitch));
    end    
    
    alpha = 1;
    
    twitch = alpha*twitch*twitch2tetanus_ratio; %*PTi(1)*0.3;
    force_temp = conv(spikeTrain,twitch);
    force = force_temp(1:length(time));
    
    force_model = Af_old(i);
    
    meanForce(i) = mean(force(3*Fs:4*Fs));
    P2PForce(i) = 1-(max(force(3*Fs:4*Fs))-min(force(3*Fs:4*Fs)))/twitch2tetanus_ratio;
    
    
end

figure(1)
plot(FR/FR_half(testingUnit),meanForce,'LineWidth',1) %; ,'Color',[0.078,0,0.831])
hold on 
plot(FR/FR_half(testingUnit),Af_old,'LineWidth',1)
xlabel('Frequency (f_{0.5})','FontSize',14)
ylabel('Force (%)','FontSize',14)



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


function spikeTrain = spikeTrainGenerator(t,Fs,freq)

spikeTrain = zeros(1,length(t));
ISI = round(1/freq*Fs);
numSpikes = round(length(t)/ISI);
index = [1:numSpikes]*ISI;
index(index>length(t)) = [];
spikeTrain(index) = 1;
spikeTrain(1) = 1;

end



