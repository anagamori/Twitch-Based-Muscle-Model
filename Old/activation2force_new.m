%--------------------------------------------------------------------------
% activationForceRelationship_test_FCR.m
% Author: Akira Nagamori
% Last update: 3/17/18
% Code descriptions
% Ojbective: fit the relationship between activation (in Hz) and force to
% the Af relationship (Af function) formulated by Loeb's model
%   - Define the minimum and peak firing rate of a given motor unit
%   - Compute f_0.5 (firing rate at which half of maximum force is achieved)
%   - Sweep firing rate from minimum (1 Hz) up to 3*f_0.5 and obtain mean
%   force level from summation of twitches
%   - Find parameters (a,b) of a scaling function that best-fits the Af
%   relationship formulated in the Loeb's model (Song et al. 2008)
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
fast_index = 89;

CT_n = 30;
CT_1 = 90;
slope_CT = (CT_1 - CT_n)/(FR_half(1)-FR_half(end));
intercept_CT = CT_1-slope_CT*FR_half(1);
CT = slope_CT*FR_half+intercept_CT;
CT = CT/1000;
RT = CT;

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

%FR = [2 5 10 14 18 22 26 30 34 38 42 46 50 54 58 62 64 68 72 76 80 84];
FR = [2 5:3:3*FR_half(testingUnit) 3*FR_half(testingUnit)];
%FR = [2 0.5*FR_half(testingUnit) 1*FR_half(testingUnit) 2*FR_half(testingUnit) 3*FR_half(testingUnit)];
% Use average twitch-tetanus ratio (0.3) for slow twitch fibers from Burke
% et al. (1974)
twitch2tetanus_ratio = 0.2;
PT = 1/twitch2tetanus_ratio;

%% Obtain non-corrected activation-force relationship
% twitch amplitude = 1
for i = 1:length(FR)
    %force = zeros(1,length(time));
    f_env = FR(i)/FR_half(testingUnit);
    spikeTrain_temp = spikeTrainGenerator(t_temp,Fs,FR(i));
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
    
    
    f_env_CT = 2*(1 - exp(-(f_env/(0.1*2.1)).^2.1));
    Af_CT = 1 - exp(-(f_env/(0.1*2.1)).^2.1);
    [twitch,T1,T2] = twitch_function(Af_old(i),Lce,CT(testingUnit),RT(testingUnit),Fs);
    
    
    %alpha = 1-exp(-((-0.5914*f_env/-1.335)^-1.335));
    if f_env < 0.5
        alpha = 1;
    else
        a1 = -0.5279;
        b1 = 0.304;
        c1 = 0.5621;
        a2 = 1.448e6;
        b2 = -74;
        c2 = 19.96;
%         a3 = 0.4441;
%         b3 = 3.611;
%         c3 = 2.904;
        %alpha = a1*exp(-((f_env-b1)/c1)^2) + a2*exp(-((f_env-b2)/c2)^2) + a3*exp(-((f_env-b3)/c3)^2);
        alpha = a1*exp(-((f_env-b1)/c1)^2) + a2*exp(-((f_env-b2)/c2)^2);
        %alpha = 0.04642*f_env^4-0.2942*f_env^3+0.5968*f_env^2-0.6327*f_env+1.201;
    end
    twitch = alpha*twitch; %*PTi(1)*0.3;
    force_temp = conv(spikeTrain,twitch);
    force = force_temp(1:length(time));
    
    force_model = Af_old(i)*PT;
    
    meanForce(i) = mean(force(3*Fs:4*Fs));
    P2PForce(i) = 1-(max(force(3*Fs:4*Fs))-min(force(3*Fs:4*Fs)))/1;
    if i == 1
        [val,loc] = min(abs((force(1*Fs+1:1.05*Fs)-meanForce(i)/2)));
        t_half(i) = (time(1*Fs+1+loc)-time(1*Fs+1));
    elseif i == 2
        [val,loc] = min(abs((force(1*Fs+1:1.05*Fs)-meanForce(i)/2)));
        t_half(i) = (time(1*Fs+1+loc)-time(1*Fs+1));
    else
        [val,loc] = min(abs((force(1*Fs+1:1.3*Fs)-meanForce(i)/2)));
        t_half(i) = (time(1*Fs+1+loc)-time(1*Fs+1));
    end
    
    figure(1)
    plot(time,force,'LineWidth',1) %,'Color',[0.078,0,0.831]
    hold on
    %    plot(time,Af_vec*PT,'LineWidth',1)
    %line([time(1) time(end)],[force_model force_model],'Color','black','LineWidth',2)
    %     hold on
    %     line([time(1) time(end)],[meanForce(i) meanForce(i)],'Color','red','LineWidth',2)
    
end

figure(2)
%plot(FR/FR_half(testingUnit),meanForce./max(meanForce)*100,'LineWidth',1,'Color',[0.078,0,0.831])
plot(FR/FR_half(testingUnit),meanForce,'LineWidth',1) %; ,'Color',[0.078,0,0.831])
hold on
plot(FR/FR_half(testingUnit),Af_old*PT,'LineWidth',1,'Color',[0.851,0.325,0.098])
xlabel('Frequency (f_{0.5})','FontSize',14)
ylabel('Force (%)','FontSize',14)


%
FR_Macefield = [5,8,10,15,20,30,50,80,100];
fusion_Macefield = [1.4,2.5,3.42,4.16,4.4,4.45,4.76,4.76,4.76]*20;
figure(3)
plot(FR/FR_half(testingUnit),P2PForce*100,'LineWidth',1,'Color',[0.078,0,0.831])
hold on 
plot(FR_Macefield/10,fusion_Macefield,'LineWidth',1)
xlabel('Frequency (Hz)','FontSize',14)
ylabel('Degree of Fusion (%)','FontSize',14)
xlim([0 4])

Error = sum(abs(meanForce./max(meanForce)*100-Af_old*100))

%% function used
function [twitch,T1,T2_temp] = twitch_function(Af,Lce,CT,RT,Fs)
a = 1.1; %1.5 for #120
b = 1.1;
T1 = CT*Lce^2+CT*Af*a;% *1.1;
T2_temp = (RT + RT*Af*b)/Lce;
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
end

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