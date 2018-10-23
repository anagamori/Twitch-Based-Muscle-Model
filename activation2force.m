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
RP_MU = 30; %range of peak tetanic tension between the smallest and largest motor units [ratio to the smallest]
b_MU = log(RP_MU)/N_MU; %coefficient to establish the exponential distribution of peak tetanic tension
PTf_MU = exp(b_MU*i_MU); %peak tetanic tension of individual motor units [ratio to the smallest]
PTi = PTf_MU./sum(PTf_MU)*F0; %correct PTf_MU to the unit of Newton [Newton]

%--------------------------------------------------------------------------
% 4) Determine the index for the largest slow-twitch fibers in a pool
F_pcsa_slow = 0.5; % fractional PSCA of slow-twitch motor units [0-1]
[~, index_slow] = min(abs(cumsum(PTi) - F0*F_pcsa_slow)); % index for the largest motor unit slow-twitch fiber in a pool

%--------------------------------------------------------------------------
% 5) Determine twithc-tetanus ratio
twitch2tetanus_ratio = 0.2; %averaged across various studies (Proske & Waite, 1974; Stephens et al. 1975; Stephens & Stuart, 1975; Burke et al. 1973; Mayer et al., 1984; Walsh et al. 1978)

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
testingUnit = 1; %index for the testing motor unit
Lce = 1; % muscle length
Y = 1;
S = 1;

FR = [2 5:3:3*FR_half(testingUnit) 3*FR_half(testingUnit)];
%FR = [2:2:3*FR_half(testingUnit)];
Af_old = length(FR);

FR_exp = [2:1:50];
Force_exp = [0.0843934948657 0.10123804947139 0.11244878344883 0.14056097931104 0.18557463705802 0.27002503920278 0.39954600637362...
    0.57977135919875 0.7318339319136 0.89515149982296 0.99087586625525 1.07532626840002 1.13724138803177 1.187888866407 ...
    1.22163488289747 1.25538089938793 1.2891269158784 1.3172391117406 1.33971748697456 1.36782968283676 1.38467423744245 ...
    1.40715261267641 1.42399716728211 1.4408417218878 1.45205877889625 1.46326951287369 1.47448024685113 1.48569098082856 ...
    1.50253553543427 1.50811244878345 1.51368936213263 1.51926627548181 1.524843188831 1.52478628155192 1.5303631949011 ...
    1.5359464312813 1.53588952400223 1.54146643735141 1.54704335070058 1.54698644342152 1.5525633567707 1.55250644949162 ...
    1.55808336284081 1.55802645556174 1.55796954828266 1.56354646163185 1.56348955435277 1.56343897010471 1.56338838585664];
FR_half_exp = 10.2;
FR_exp = FR_exp/FR_half_exp;
Force_exp = Force_exp./(Force_exp(end));
% FR_exp = [2 5 8 10 15 20 30 50 80 100];
% Force_exp = [4.44 13.5 30.52 46.26 64.34 76.32 87.24 92.07 92.26 94.53];
% Force_exp = Force_exp/Force_exp(end);
% FR_half_exp = 10.2;
% FR_exp_int = 2:0.1:100;
% Force_exp_int = spline(FR_exp,Force_exp,FR_exp_int);
% FR_exp_int = FR_exp_int/FR_half_exp;

%plot(FR_exp_int,Force_exp_int)
%%
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
    
    %----------------------------------------------------------------------
    % Compute the desired frequency-force relationship
    [~, index_FR_half] = min(abs(FR_exp - f_env));
    Force_desired = Force_exp(index_FR_half);
    %----------------------------------------------------------------------
    % Compute a time-series of Af
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
            Af_temp = Af_slow_function(f_eff,Lce,Y);
        else
            Af_temp = Af_fast_function(f_eff,Lce,S);
        end
        Af_vec(t) = Af_temp;
    end
    %----------------------------------------------------------------------
    % Create a tiwthc profile
    T1 = CT(testingUnit)+CT(testingUnit)*(f_env/2)^2;
    T2_temp = RT(testingUnit)+RT(testingUnit)*(f_env/4)^2;
    
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
    alpha = 0.2;
    twitch = alpha*twitch; %*PTi(1)*0.3;
    Force_temp = conv(spikeTrain,twitch);
    Force = Force_temp(1:length(time));
    %Force = 1./(1+exp(-5*(Force-1)));
    Force = Force.^2./(Force.^2+0.5^2);
    %Force = 0.5*(1+tanh((Force-0.5)/1));
    mean_Force = mean(Force(2*Fs:3*Fs));
    %----------------------------------------------------------------------
    % Find a coefficient, alpha, to fit model force to the desired force
%     alpha = 1*rand(1);
%     for j = 1:100
%         Force_test = alpha*mean_Force;
%         if Force_test < Force_desired
%             alpha = alpha + 0.5./2.^(j-2);
%         else
%             alpha = alpha - 0.5./2.^(j-2);
%         end
%     end
    
    figure(1)
    plot(time,Force)
    hold on 
    plot([time(1) time(end)],[Force_desired Force_desired])
    
    Force_model = Af_old(i);
    
    meanForce(i) = mean(Force(3*Fs:4*Fs));
    P2PForce(i) = max(Force(3*Fs:4*Fs))-min(Force(3*Fs:4*Fs));
    
    alpha_vec(i) = alpha;
    
end

freq = FR/FR_half(testingUnit);
P2PForce = (1-P2PForce/P2PForce(1))*100;

figure(2)
plot(FR/FR_half(testingUnit),meanForce,'LineWidth',1) %; ,'Color',[0.078,0,0.831])
hold on
plot(FR/FR_half(testingUnit),Af_old,'LineWidth',1)
hold on 
plot(FR_exp,Force_exp,'LineWidth',1)
xlabel('Frequency (f_{0.5})','FontSize',14)
ylabel('Force (%)','FontSize',14)
xlim([0 3.1])

figure(3)
plot(FR/FR_half(testingUnit),P2PForce,'LineWidth',1) %; ,'Color',[0.078,0,0.831])
hold on 

figure(4)
plot(FR/FR_half(testingUnit),meanForce,'LineWidth',1) %; ,'Color',[0.078,0,0.831])
hold on

% figure(4)
% plot(FR/FR_half(testingUnit),alpha_vec,'LineWidth',1)

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



