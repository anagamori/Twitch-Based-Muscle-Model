%--------------------------------------------------------------------------
% activation2force_new_fit_v2.m
% Author: Akira Nagamori
% Last update: 7/11/18
% Code descriptions
% Ojbective:
%--------------------------------------------------------------------------

close all
clear all
clc
%--------------------------------------------------------------------------
% motor unit parameters
N_MU = 120; % number of motor units
i_MU = 1:N_MU; % index for motor units

%--------------------------------------------------------------------------
L0 = 3;
density = 1.06; % muscle density [g/cm^3]
mass = 0.0055; % muscle mass [kg]
PCSA = (mass*1000)/(density*L0); % PCSA of muscle
sigma = 22.4; % specific tension
F0 = PCSA*sigma;

%--------------------------------------------------------------------------
% Peak tetanus parameter
RP_MU = 25; %range of peak tension across motor untis in unit of fold
b_MU = log(RP_MU)/N_MU; %coefficient to establish a range of twich force values
P_MU = exp(b_MU*i_MU); %force generated by a motor unit as a function of its recruitment threshold

PTi = P_MU./sum(P_MU)*F0;
a_twitch = 0.014061531587007;
b_twitch = 0.030319762726547;
Pi_MU = a_twitch*exp(b_twitch*i_MU);

F_pcsa_slow = 0.5; % fractional PSCA of slow-twitch motor units (0-1)
[~, index_slow] = min(abs(cumsum(PTi) - F0*F_pcsa_slow)); % index for the largest motor unit clacified as slow-twitch
%--------------------------------------------------------------------------
% Recruitment threshold of each unit
Ur = 0.6; % recruitment threshold for the lastly recruited motor unit
Ur_1 = 0.01; % reruitment threshold for the first unit
b_Ur = log(Ur/Ur_1)/N_MU;
U_th = exp(b_Ur*i_MU)/100;
%--------------------------------------------------------------------------
% Minimum and maximum firing rates of each unit
MFR1_MU = 8; %minimum firing rate of first unit
MFRn_MU = 14; %minimum firing rate of last unit
RTEn_MU = U_th(end)-Ur_1;  %recruitment threshold of last unit
MFR_MU = MFR1_MU + (MFRn_MU-MFR1_MU) * ((U_th-Ur_1)./RTEn_MU);
PFR_MU = 4*MFR_MU; %peak firing rate
FR_half = PFR_MU./2; % firing rate at which half of maximum tension is achieved

%--------------------------------------------------------------------------
% Contraction time and half-relaxation time
CT_n = 30;
CT_1 = 90;
slope_CT = (CT_1 - CT_n)/(FR_half(1)-FR_half(end));
intercept_CT = CT_1-slope_CT*FR_half(1);
CT = slope_CT*FR_half+intercept_CT;
CT = CT/1000;
RT = CT;

%--------------------------------------------------------------------------
Fs = 1000;
time = 0:1/Fs:5;
t_temp = 0:1/Fs:3;
%--------------------------------------------------------------------------
% simulation parameters
Lce = 1;
Y = 1;
S = 0.96;

%--------------------------------------------------------------------------
% Find an index of representative slow twitch fiber
% From Burke et al. (1974), fibers types between slow and fast can be differentiable based on the length of contraction time
% A unit that has contraction time closest to mean contraction time of all
% slow twitch fibers is defined as a representative unit of slow twitch.
mean_CT_slow = mean(CT(1:index_slow));
[~, index_slow_rep] = min(abs(CT(1:index_slow) - mean_CT_slow));

testingUnit =  120; %index_slow_rep;
%FR_test = [12 16 20 24 28 32 36 40 44 48];
FR_test = [20 24 28 32 36 40 44 48 52 56 60 64 68 72 76 80 84];
for j = 1:length(FR_test)
    FR = FR_test(j); %round(3*FR_half(testingUnit)); %round(3*FR_half(testingUnit)) %round(2.7*FR_half(testingUnit)); % [2 5 10 15 20 25 30 35 40 45 48];
    f_env = FR/FR_half(testingUnit);
    
    % Use average twitch-tetanus ratio (0.2) for slow twitch fibers from
    % Macefield et al. (1993)
    twitch2tetanus_ratio = 0.2;
    PT = 1/twitch2tetanus_ratio;
    %--------------------------------------------------------------------------
    % Generate spike train
    spikeTrain_temp = spikeTrainGenerator(t_temp,Fs,FR);
    spikeTrain = [zeros(1,1*Fs) spikeTrain_temp zeros(1,1*Fs)];
    
    %% Obtain reference activation-force relationship
    %--------------------------------------------------------------------------
    % Parameter initialization
    f_env_2 = [zeros(1,1*Fs) f_env*ones(1,length(spikeTrain_temp)) zeros(1,1*Fs)];
    f_int = 0;
    f_eff = 0;
    f_eff_dot = 0;
    f_eff_vec = zeros(1,length(f_env_2));
    Af_vec = zeros(1,length(f_env_2));
    
    %--------------------------------------------------------------------------
    % Activation-frequency relationship from Song et al. 2008
    if testingUnit <= index_slow
        Af_old = Af_slow_function(f_env,Lce,Y);
    else
        Af_old = Af_fast_function(f_env,Lce,S);
    end
    
    %--------------------------------------------------------------------------
    % Reference rise and fall time of muscle force from Song et al. (2008)
    for t = 1:length(f_env_2)
        
        if f_eff_dot >= 0
            if testingUnit <= index_slow
                T_f = 0.0343*Lce^2 + 0.0227*f_env;
            else
                T_f = 0.0206*Lce^2 + 0.0136*f_env;
            end
        else
            if testingUnit <= index_slow
                T_f = (0.047+0.0252*Af_old)/Lce;
            else
                T_f = (0.0282+0.0151*Af_old)/Lce;
            end
        end
        f_int_dot = (f_env_2(t) - f_int)/T_f;
        f_int = f_int_dot*1/Fs + f_int;      
        f_eff_dot = (f_int - f_eff)/T_f;
        f_eff = f_eff_dot*1/Fs + f_eff;
        f_eff_vec(t) = f_eff;
        if testingUnit <= index_slow
            Af_vec(t) = Af_slow_function(f_eff,Lce,Y);
        else
            Af_vec(t) = Af_fast_function(f_eff,Lce,S);
        end
    end
    %--------------------------------------------------------------------------
    % Calculate rise (t_0-50) and fall time (t_100-50)
    % Rise time: time from stimulation to 50% of peak force (Brown & Loeb, 2000, IV)
    % Fall time: time from end of stimulation to 50% of peak of the force that was present at the end of stimulation (Brown & Loeb, 2000, IV)
    maxForce_1 = max(Af_vec);
    [~,loc_rise_1] = min(abs((Af_vec(1*Fs+1:1.3*Fs)-maxForce_1/2)));
    t_rise_time_1 = (time(1*Fs+1+loc_rise_1)-time(1*Fs+1));
    [~,loc_fall_1] = min(abs((Af_vec(4*Fs+1:4.3*Fs)-maxForce_1/2)));
    t_fall_time_1 = (time(4*Fs+1+loc_fall_1)-time(4*Fs+1));
    
    mean_Force_1 = mean(Af_vec(2*Fs:3*Fs)*PT);
    %--------------------------------------------------------------------------
    %% Obtain new activation-force relationship
    %--------------------------------------------------------------------------
    % Create twitch profile
    a = 1; %2*rand(1);
    b = 1; %2*rand(1);
    
    T1 = CT(testingUnit)*Lce^2 + CT(testingUnit)*Af_old*a;
    T2_temp = (RT(testingUnit)+ RT(testingUnit)*Af_old*b)/Lce; %;
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
    
    %--------------------------------------------------------------------------
    % Compute resulting force with a given twithc profile
    alpha = 2*rand(1);
    for i = 1:100
        Force_temp = conv(spikeTrain,twitch*alpha);
        Force = Force_temp(1:length(time));
        mean_Force_2 = mean(Force(2*Fs:3*Fs));
        if mean_Force_2 < mean_Force_1
            alpha = alpha + 0.6./2.^(i-2);
        else
            alpha = alpha - 0.6./2.^(i-2);
        end
    end
    %--------------------------------------------------------------------------
    % Calculate rise (t_0-50) and fall time (t_100-50)
    maxForce_2 = max(Force);
    spike_time = find(spikeTrain == 1);
    end_stimulation = spike_time(end);
    endForce_2 = Force(end_stimulation);
    [~,loc_rise_2] = min(abs((Force(1*Fs+1:1.3*Fs)-maxForce_2/2)));
    t_rise_time_2 = (time(1*Fs+1+loc_rise_2)-time(1*Fs+1));
    [~,loc_fall_2] = min(abs((Force(end_stimulation:end_stimulation+0.3*Fs)-endForce_2/2)));
    t_fall_time_2 = (time(end_stimulation+loc_fall_2)-time(end_stimulation ));
    
    maxForce_2_vec(j) = maxForce_2;
    %--------------------------------------------------------------------------
    % Calculate degree of fusion (Macefield et al. (1993))
    p2p_Force = max(Force(2*Fs:3*Fs))-min(Force(2*Fs:3*Fs));
    fusion(j) = (1 - p2p_Force/1)*100;
    %--------------------------------------------------------------------------
    % Plot results
    figure(1)
    plot(time,Af_vec*PT,'LineWidth',1)
    hold on
    plot(time,Force,'LineWidth',1)
    
    % display results
    
    t_rise_time_1_vec(j) = t_rise_time_1;
    t_rise_time_2_vec(j) = t_rise_time_2;
    t_fall_time_1_vec(j) = t_fall_time_1;
    t_fall_time_2_vec(j) = t_fall_time_2;
    
    alpha_vec(j) = alpha;
    
end

% Reference frequency-force relationship from Bellemare 1983
% FR_half_reference = 14.4;
% FR_reference = [3 5 8 10 15 20 30 40 50 60 80];
% Force_reference = [6.34 7.3 13.7 22.7 55.1 71.5 87.2, 93.7 96 98.9 100];

figure(2)
plot(FR_test,fusion,'LineWidth',1)
xlabel('Frequency (Hz)')
ylabel('Fusion')

x = FR_test/FR_half(testingUnit);
x = x';
%x = [0.5;x];
y = alpha_vec;
y = y';
%y = [0.85; y];
%% function used
function Af = Af_slow_function(f_eff,L,Y)
a_f = 0.56;
n_f0 = 1.9;
n_f1 = 5;
n_f = n_f0 + n_f1*(1/L-1);
Af = 1 - exp(-(Y*f_eff/(a_f*n_f))^n_f);
end

function Af = Af_fast_function(f_eff,L,S)
a_f = 0.56;
n_f0 = 1.9;
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