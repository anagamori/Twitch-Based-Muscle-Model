%--------------------------------------------------------------------------
% test_v3.m
% Author: Akira Nagamori
% Last update: 4/18/18
% Code descriptions
% Ojbective:
%   Test activation-force and activation-fusion relationships
%--------------------------------------------------------------------------
close all
clear all
clc

%--------------------------------------------------------------------------
N_MU = 120; % number of motor units
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
FR_half = PFR_MU./2;

CT_n = 20;
FR_half_n = FR_half(end);
CT = 1.5*(CT_n*FR_half_n)./FR_half;
CT = CT - (CT(end)-CT_n);
CT = CT/1000;
RT = CT+0.03;

%--------------------------------------------------------------------------
a = 0.56;
n = 2.1;
%--------------------------------------------------------------------------
Fs = 1000;
time = 0:1/Fs:4;
time_stim = 0:1/Fs:2;
time_twitch = 0:1/Fs:3;

testingUnit = 120;
FR = [0:0.01:3]*FR_half(testingUnit);
p2p_force = zeros(1,length(FR));

for i = 1:length(FR)
    test_freq = FR(i);
    f_env = test_freq/FR_half(testingUnit);
    spikeTrain_temp = spikeTrainGenerator(time_stim,Fs,test_freq);
    spikeTrain = [zeros(1,Fs) spikeTrain_temp zeros(1,Fs)];
    
    alpha_CT = 1-exp(-(f_env/a)^n);
    T_c = CT(testingUnit)+CT(testingUnit)*alpha_CT;
    T_hr = RT(testingUnit)+RT(testingUnit)*alpha_CT;
    
    k = log(2)/(T_hr-T_c-T_c*log(T_hr/T_c));
    m = k*T_c;
    p = exp(-k*T_c*(log(T_c)-1));
    twitch = p*time_twitch.^m.*exp(-k*time_twitch);
    
    force_temp = conv(spikeTrain,twitch);
    force = force_temp(1:length(time));
    
    p2p_force(i) = max(force(2*Fs:3*Fs))-min(force(2*Fs:3*Fs));
    figure(1)
    plot(time,force)
    hold on
end

figure(2)
plot(FR./FR_half(testingUnit),p2p_force)


%%
FR = 0:0.01:3;
beta = zeros(1,length(FR));

for j = 1:length(FR)
    if FR(j) <= 0.3
        beta(j) = 1;
    elseif FR(j) > 0.3 && FR(j) <= 1
        beta(j) = -exp(-6*(FR(j)-0.3))*0.5+1.5;
    else
        beta(j) = (1.5-0.3)*exp(-3*(FR(j)-1))+0.3;
    end
end

figure(3)
plot(FR,beta)
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
