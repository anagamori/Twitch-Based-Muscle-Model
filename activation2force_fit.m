%==========================================================================
% activation2force_fit.m
% Author: Akira Nagamori
% Last update: 10/25/18
% Descriptions: find parameters, a and b, for non-linear scaling of force
% output, for each motor unit
%==========================================================================
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
% 5) Determine the recruitment threshold of motor units [0 1]
Ur = 0.6; % recruitment threshold for the last-recruited motor unit (i.e., the largest motor unit)
Ur_1 = 0.01; % reruitment threshold for the first unit
b_Ur = log(Ur/Ur_1)/N_MU; %coefficient to establish the exponential distribution of recruitment threshold
U_th = exp(b_Ur*i_MU)/100; % recruitment threshold for each motor unit [0 1]

%--------------------------------------------------------------------------
% 6) Determine the contraction time and half-reflaxation time
% Use the formulation given by Fuglevand et al. (1993)
CT_n = 30; %contraction time of first unit [ms]
CT_1 = 90; %contraction time of lasdt unit [ms]
RC_MU = CT_1/CT_n; %range of contraction time
c = log(RP_MU)/log(RC_MU); %scaling factor 
CT = CT_1*(1./PTf_MU).^(1/c); %assign contraction 
CT = CT/1000; %contraction time [s]
RT = CT; %half-relaxation time [s]

%--------------------------------------------------------------------------
% 7) Determine FR_half (frequency at which half the maximum tetanic tension is achieved)
% Use the relationship between contraction time and FR_half drived by Kernell et al. (1983)
FR_half = zeros(1,N_MU);
for n = 1:N_MU
    if n <= index_slow
        FR_half(n) = 1/(CT(n)*1.5); %slow-twitch fibers
    else
        FR_half(n) = 1/(CT(n)*1.8); %fast-twitch fibers
    end
end
%--------------------------------------------------------------------------
% 8) Determine minimum and peak firing rate of each motor unit 
% Use the formulation given in Song et al. (2008)
% Minimum firing rate = half the FR_half
% Peak firing rate = twitch the FR_half
MFR_MU = FR_half./2;
PFR_MU = FR_half*2;

%==========================================================================
% Parameters for sag (Brown & Loeb, 2000 IV)
T_S = 0.043;
a_S1 = 1.76;
a_S2 = 0.96;
S = 0;
S_temp = 0.96;

%--------------------------------------------------------------------------
% Parameters for yielding (Brown & Loeb, 2000 IV)
Y = 1;

%==========================================================================
% Frequency to force relationship from Rack & Westbury 1969
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

%==========================================================================
% Simulation
Fs = 1000; %sampling frequency
time = 0:1/Fs:5; %simulation time
time_active = 0:1/Fs:3; %time window in which a motor unit discharges

testingUnit = 300; %index for the testing motor unit
Lce = 1; % muscle length

FR = [2 5:3:3*FR_half(testingUnit) 3*FR_half(testingUnit)]; % testing firing rates

S_vec = zeros(1,length(time)); % vector for sag parameter, S
Af = length(FR);
meanForce = length(FR); %mean force 
P2PForce = length(FR);
    
%==========================================================================
% Start simulation
for i = 1:length(FR)        
    %----------------------------------------------------------------------
    % Compute the value of Af (activation-frequency relationship) from
    % Tsianos et al. (2012)
    f_env = FR(i)/FR_half(testingUnit); % normalized frequency Brown & Loeb (1999)
    if testingUnit <= index_slow
        Af(i) = Af_slow_function(f_env,Lce,Y); %activation-frequency relationship for slow-twitch fibers    
    else
        Af(i) = Af_fast_function(f_env,Lce,S_temp); %activation-frequency relationship for fast-twitch fibers          
    end
    
    %----------------------------------------------------------------------
    % Generate a spike train at the given frequency 
    spikeTrain_temp = spikeTrainGenerator(time_active,Fs,FR(i));
    spikeTrain = [zeros(1,1*Fs) spikeTrain_temp zeros(1,1*Fs)];
    
    %----------------------------------------------------------------------
    % Initialize vectors for the follwoing sections
    f_env_2 = [zeros(1,1*Fs) f_env*ones(1,length(spikeTrain_temp)) zeros(1,1*Fs)];
    f_int = 0;
    f_eff = 0;
    f_eff_dot = 0;
    f_eff_vec = zeros(1,length(f_env_2));
    Af_vec = zeros(1,length(f_env_2));
    
    %----------------------------------------------------------------------
    % Compute a time-series of Af values from Tsianos et al. (2012)
    for t = 1:length(f_env_2)
        if f_eff_dot >= 0
            if testingUnit <= index_slow
                T_f = 0.0343 + 0.0227*f_env;
            else
                T_f = 0.0206 + 0.0136*f_env;
            end
        else
            if testingUnit <= index_slow
                T_f = (0.047+0.0252*Af(i))/Lce;
            else
                T_f = (0.0282+0.0151*Af(i))/Lce;
            end
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
    
    %======================================================================
    % Testing a new model
    %----------------------------------------------------------------------
    % Create a tiwthc profile
    T1 = CT(testingUnit)+CT(testingUnit)/2*f_env; % contraction time
    T2_temp = RT(testingUnit)+RT(testingUnit)/4*f_env; % half-relaxation time
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
    % Determine twithc-tetanus ratio
    if testingUnit <= index_slow
        twitch2tetanus_ratio = 0.2; %slow-twitch fibers
    else
        twitch2tetanus_ratio = 0.4; %fast-twitch fibers
    end
    
    % 0.2 for average twitch-tetanus ratio for slow tiwtich
    % 0.457 for average tiwtch-tetnaus ratio for fast twitch (both FF and FR):
    % (Stephens & Stuart, 1975b; Burke et al. 1973; Devanandan et al.,
    % 1965; Dum & Kennedy 1980; Mayer et al., 1984; Walsh et al., 1978)
    twitch = twitch2tetanus_ratio*twitch;
    
    %----------------------------------------------------------------------
    % Generate force output by convolving the twitch profile with the spike
    % train
    Force_temp = conv(spikeTrain,twitch);
    Force = Force_temp(1:length(time));
    
    %----------------------------------------------------------------------
    % Apply sag if the unit is fast-twitch
    if testingUnit > index_slow
        for t = 1:length(time)
            if Force(t) < 0.1
                a_S = a_S1;
            else
                a_S = a_S2;
            end
            S_dot = (a_S - S)/T_S;
            S = S_dot*1/Fs + S;
            
            S_vec(t) = S;
        end
        Force = Force.*S_vec;
    end
    
    %----------------------------------------------------------------------
    % Non-linear scaling of output force
    a = 2.2; % unit 1: 2.2
    b = 0.48; % unit 1:0.3
    Force = Force.^a./(Force.^a+b^2);    
          
    %----------------------------------------------------------------------
    % Plot output
    figure(1)
    plot(time,Force) % new model
    hold on
    plot(time,Af_vec) % Tsianos et al. (2012)
      
    %----------------------------------------------------------------------
    % Calculate output variables
    meanForce(i) = mean(Force(3*Fs:4*Fs)); %mean force 
    P2PForce(i) = max(Force(3*Fs:4*Fs))-min(Force(3*Fs:4*Fs)); %peak-to-peak amplitude in force   
    
end

%--------------------------------------------------------------------------
P2PForce = (1-P2PForce/P2PForce(1))*100; %convet peak-to-peak amplitude into fusion index (Macefield et al., 1993)

%----------------------------------------------------------------------
% Plot frequency-force relationship 
figure(2)
plot(FR/FR_half(testingUnit),meanForce,'LineWidth',1) %; ,'Color',[0.078,0,0.831])
hold on
plot(FR/FR_half(testingUnit),Af,'LineWidth',1)
hold on
plot(FR_exp,Force_exp,'LineWidth',1)
xlabel('Frequency (f_{0.5})','FontSize',14)
ylabel('Force (%)','FontSize',14)
xlim([0 3.1])
hold on
plot([0 3.1],[0.5 0.5],'k','LineWidth',1)
plot([1 1],[0 1],'k','LineWidth',1)
legend('New','Song','Rack')

%----------------------------------------------------------------------
% Plot fusion index
figure(3)
plot(FR/FR_half(testingUnit),P2PForce,'LineWidth',1) %; ,'Color',[0.078,0,0.831])
hold on

%----------------------------------------------------------------------
% Plot force vs fusion (Macefield et al., 19993)
figure(4)
plot(meanForce*100,P2PForce,'LineWidth',1)
hold on
plot(1:100,1:100,'--','LineWidth',1)

%==========================================================================
% Function used 
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



