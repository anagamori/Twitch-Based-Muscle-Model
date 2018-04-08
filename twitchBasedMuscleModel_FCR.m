%--------------------------------------------------------------------------
% twitchBasedMuscleModel_FCR.m
% Author: Akira Nagamori
% Last update: 4/5/17
% ---Code descriptions-----
% Ojbective: find twitch amplitudes for each motor unit
%
%--------------------------------------------------------------------------
function output = twitchBasedMuscleModel_FCR(amplitude)

%--------------------------------------------------------------------------
alpha = 3.1*pi/180;
L0 = 5.1; % optimal muscle length [cm]
Lse_slack = 27.1;
L0T = Lse_slack*1.05;
Lce_initial = 5.1;
Lse_initial = 27.1;
Lmt = Lce_initial+cos(alpha)*Lse_initial;
[Lce,Lse,Lmax] =  InitialLength_function(L0,alpha,Lse_slack,Lce_initial,Lse_initial);

%--------------------------------------------------------------------------
% F0 parameters
density = 1.06; % muscle density [g/cm^3]
mass = 0.02; % muscle mass [kg]
PCSA = (mass*1000)/(density*L0); % PCSA of muscle
sigma = 31.8; % specific tension
F0 = PCSA*sigma;
N_MU = 300; % number of motor units
i_MU = 1:N_MU; % index for motor units

%--------------------------------------------------------------------------
% Peak tetanus parameter
RP_MU = 25; %range of peak tension across motor untis in unit of fold
b_MU = log(RP_MU)/N_MU; %coefficient to establish a range of twich force values
P_MU = exp(b_MU*i_MU); %force generated by a motor unit as a function of its recruitment threshold

PTi = P_MU./sum(P_MU)*F0;
a_twitch = 0.016860913109322;
b_twitch = 0.012195538822182;
Pi_MU = a_twitch*exp(b_twitch*i_MU);
%--------------------------------------------------------------------------
% Recruitment parameters
%   Find recruitment threshold for individual units using exponential fit
F_pcsa_slow = 0.5; % fractional PSCA of slow-twitch motor units (0-1)
[~, index_slow] = min(abs(cumsum(PTi) - F0*F_pcsa_slow)); % index for the largest motor unit clacified as slow-twitch
Ur = 0.6; % recruitment threshold for the last recruited motor unit
Ur_1 = 0.01; % reruitment threshold for the first unit
f_RT = fit([1 N_MU]',[Ur_1 Ur]','exp1');
coeffs_f_RT = coeffvalues(f_RT);
U_th = coeffs_f_RT(1)*exp(coeffs_f_RT(2)*i_MU); % the resulting recruitment threshold for individual units
index_fast = 239;
%--------------------------------------------------------------------------
% Firing rate parameters
%   Find peak firing rates of individual units
MFR1_MU = 8; %minimum firing rate of first unit
MFRn_MU = 14; %minimum firing rate of last unit
RTEn_MU = U_th(end)-Ur_1;  %recruitment threshold of last unit
MFR_MU = MFR1_MU + (MFRn_MU-MFR1_MU) * ((U_th-Ur_1)./RTEn_MU);
PFR_MU = 4*MFR_MU; %peak firing rate
FR_half = PFR_MU./2; % firing rate at which half of maximum tension is achieved
g_e = (PFR_MU(end) - MFR_MU(end))/(1-U_th(end));
%--------------------------------------------------------------------------
% Contraction time and half relaxation time
%    Find contraction time and half relaxation time of individual units
%   Using Loeb's formulation that contraction time of individual units is
%   proportional to 1/FR_half (Brown & Loeb 2000 IV; Botterman et al., 1996)
CT_n = 20;
FR_half_n = FR_half(end);
CT = 1.5*(CT_n*FR_half_n)./FR_half;
CT = CT - (CT(end)-CT_n);
CT = CT/1000;
RT = CT;

%--------------------------------------------------------------------------
cv_MU = 0.1; %ISI variability as per coefficient of variation (=mean/SD)

%--------------------------------------------------------------------------

Vce = 0;
Y = 0;
U_eff = 0;
Force = 0;
Fs = 40000;
h = 1/Fs;
time = 0:1/Fs:10;

amp = amplitude;
U = [zeros(1,1*Fs) (amp/2)*(0:1/Fs:2) amp*ones(1,length(time)-3*Fs-1)];


FR = zeros(1,N_MU);
FR_mat = zeros(N_MU,length(time));
f_env = zeros(1,N_MU);
S = zeros(1,N_MU);
spike_train = zeros(N_MU,length(time));
spike_time = zeros(1,N_MU);
spike_time_previous = zeros(1,N_MU);
force = zeros(N_MU,length(time));


MuscleVelocity = zeros(1,length(time));
MuscleLength = zeros(1,length(time));
MuscleLength(1) = Lce*L0/100;

Force_vec = zeros(1,length(time));
Tendon_Force_vec = zeros(1,length(time));
Vce_vec = zeros(1,length(time));
Lce_vec = zeros(1,length(time));
Lse_vec = zeros(1,length(time));
U_eff_vec = zeros(1,length(time));

for t = 1:length(time)
    if U(t) >= U_eff
        T_U = 0.03;
    else
        T_U = 0.15;
    end
    U_eff_dot = (U(t)-U_eff)/T_U;
    U_eff = U_eff_dot*1/Fs+U_eff;
    
    %----------------------------------------------------------------------
    % Yield
    Y = yield_function(Y,Vce,Fs);
    
    %----------------------------------------------------------------------
    for n = 1:N_MU
        FR(n) = g_e.*(U_eff - U_th(n)) + MFR_MU(n);
        if FR(n) < MFR_MU(n)
            FR(n) = 0;
        elseif FR(n) > PFR_MU(n)
            FR(n) = PFR_MU(n);
        end
        FR_mat(n,t) = FR(n);
        f_env(n) = FR(n)/FR_half(n);
        if n > index_fast
            S(n) = sag_function(S(n),f_env(n),Fs);
        end
    end
    
    %----------------------------------------------------------------------
    if t > 1
        index_1 = i_MU(FR >= MFR_MU & FR_mat(:,t-1)' == 0);
        index_2 = i_MU(FR >= MFR_MU & spike_time == t);
        index = [index_1 index_2];
        
        for j = 1:length(index)
            spike_train_temp = zeros(1,length(time));
            index_MU = index(j);
            noise_FR = FR(index_MU);
            if ~any(spike_train(index_MU,:)) % initial time
                spike_train(index_MU,t) = 1;
                spike_train_temp(t) = 1;
                mu = 1/FR(index_MU);
                Z = randn(1);
                Z(Z>3.9) = 3.9;
                Z(Z<-3.9) = -3.9;
                noise = 1/noise_FR*cv_MU*Z;
                spike_time_temp = (mu + noise)*Fs;
                if spike_time_temp < 2*1000/Fs
                    spike_time_temp = 2;
                end
                spike_time(index_MU) = round(spike_time_temp) + t;
                
                if n <= index_fast
                    Af = Af_slow_function(f_env(index_MU),Lce,Y);
                    Af_cor = Af_slow_correction_function(f_env(index_MU),Lce,Y);
                else
                    Af = Af_fast_function(f_env(index_MU),Lce,S(index_MU));
                    Af_cor = Af_fast_correction_function(f_env(index_MU),Lce,S(index_MU));
                end
                
                [twitch,~,~] = twitch_function(Af,Lce,CT(index_MU),RT(index_MU),Fs);
                force_temp = conv(spike_train_temp,Pi_MU(index_MU)*twitch*Af_cor);
                force(index_MU,:) = force(index_MU,:) + force_temp(1:length(time));
                spike_time_previous(index_MU) = t;
            else
                if spike_time(index_MU) == t
                    spike_train(index_MU,t) = 1;
                    spike_train_temp(t) = 1;
                    mu = 1/FR(index_MU);
                    Z = randn(1);
                    Z(Z>3.9) = 3.9;
                    Z(Z<-3.9) = -3.9;
                    noise = 1/noise_FR*cv_MU*Z;
                    spike_time_temp = (mu + noise)*Fs;
                    if spike_time_temp < 2*1000/Fs
                        spike_time_temp = 2;
                    end
                    spike_time(index_MU) = round(spike_time_temp) + t;
                    
                    ISI = (t-spike_time_previous(index_MU))/Fs;
                    FR_temp = 1/ISI;
                    if index_MU <= index_fast
                        Af = Af_slow_function(FR_temp/FR_half(index_MU),Lce,Y);
                        Af_cor = Af_slow_correction_function(FR_temp/FR_half(index_MU),Lce,Y);
                    else
                        Af = Af_fast_function(FR_temp/FR_half(index_MU),Lce,S(index_MU));
                        Af_cor = Af_fast_correction_function(FR_temp/FR_half(index_MU),Lce,S(index_MU));
                    end                   
                    [twitch,~,~] = twitch_function(Af,Lce,CT(index_MU),RT(index_MU),Fs);
                    force_temp = conv(spike_train_temp,Pi_MU(index_MU)*twitch*Af_cor);
                    force(index_MU,:) = force(index_MU,:) + force_temp(1:length(time));
                    spike_time_previous(index_MU) = t;
                elseif FR_mat(index_MU,t-1) == 0
                    spike_train(index_MU,t) = 1;
                    spike_train_temp(t) = 1;
                    
                    mu = 1/FR(index_MU);
                    Z = randn(1);
                    Z(Z>3.9) = 3.9;
                    Z(Z<-3.9) = -3.9;
                    noise = 1/noise_FR*cv_MU*Z;
                    spike_time_temp = (mu + noise)*Fs;
                    if spike_time_temp < 2*1000/Fs
                        spike_time_temp = 2;
                    end
                    spike_time(index_MU) = round(spike_time_temp) + t;
                    if n <= index_fast
                        Af = Af_slow_function(f_env(index_MU),Lce,Y);
                        Af_cor = Af_slow_correction_function(f_env(index_MU),Lce,Y);
                    else
                        Af = Af_fast_function(f_env(index_MU),Lce,S(index_MU));
                        Af_cor = Af_fast_correction_function(f_env(index_MU),Lce,S(index_MU));
                    end
                    [twitch,~,~] = twitch_function(Af,Lce,CT(index_MU),RT(index_MU),Fs);
                    force_temp = conv(spike_train_temp,Pi_MU(index_MU)*twitch*Af_cor);
                    force(index_MU,:) = force(index_MU,:) + force_temp(1:length(time));
                    
                    spike_time_previous(index_MU) = t;
                end
            end
        end
        
        FL_slow = FL_slow_function(Lce);
        FL_fast = FL_fast_function(Lce);
        if Vce > 0
            FV_slow = FV_ecc_slow_function(Lce,Vce);
            FV_fast = FV_ecc_fast_function(Lce,Vce);
        else
            FV_slow = FV_con_slow_function(Lce,Vce);
            FV_fast = FV_con_fast_function(Lce,Vce);
        end
        
        force(1:index_slow,t) = force(1:index_slow,t)*FL_slow*FV_slow;
        force(index_slow+1:end,t) = force(index_slow+1:end,t)*FL_fast*FV_fast;
        
        FP1 = F_pe_1_function(Lce/Lmax,Vce);
        FP2 = F_pe_2_function(Lce);
        if FP2 > 0
            FP2 = 0;
        end
        
        Force = sum(force(:,t)) + FP1*F0 + FP2*F0;
    end
    if Force < 0
        Force = 0;
    end
    ForceSE = F_se_function(Lse) * F0;
    
    k_0 = h*MuscleVelocity(t);
    l_0 = h*((ForceSE*cos(alpha) - Force*(cos(alpha)).^2)/(mass) ...
        + (MuscleVelocity(t)).^2*tan(alpha).^2/(MuscleLength(t)));
    k_1 = h*(MuscleVelocity(t)+l_0/2);
    l_1 = h*((ForceSE*cos(alpha) - Force*(cos(alpha)).^2)/(mass) ...
        + (MuscleVelocity(t)+l_0/2).^2*tan(alpha).^2/(MuscleLength(t)+k_0/2));
    k_2 = h*(MuscleVelocity(t)+l_1/2);
    l_2 = h*((ForceSE*cos(alpha) - Force*(cos(alpha)).^2)/(mass) ...
        + (MuscleVelocity(t)+l_1/2).^2*tan(alpha).^2/(MuscleLength(t)+k_1/2));
    k_3 = h*(MuscleVelocity(t)+l_2);
    l_3 = h*((ForceSE*cos(alpha) - Force*(cos(alpha)).^2)/(mass) ...
        + (MuscleVelocity(t)+l_2).^2*tan(alpha).^2/(MuscleLength(t)+k_2));
    MuscleLength(t+1) = MuscleLength(t) + 1/6*(k_0+2*k_1+2*k_2+k_3);
    MuscleVelocity(t+1) = MuscleVelocity(t) + 1/6*(l_0+2*l_1+2*l_2+l_3);
    
    % normalize each variable to optimal muscle length or tendon length
    Vce = MuscleVelocity(t+1)/(L0/100);
    Lce = MuscleLength(t+1)/(L0/100);
    Lse = (Lmt - Lce*L0*cos(alpha))/L0T;
    
    U_eff_vec(t) = U_eff;
    
    Force_vec(t) = Force; % muscle force
    Tendon_Force_vec(t) = ForceSE; % tendon force
    Lse_vec(t) = Lse; % normalized tendon length
    Lce_vec(t) = Lce; % normalized muscle length
    Vce_vec(t) = Vce; % normalized muscle excursion velocity
    
end

output.Force = force;
output.Tendon_Force = Tendon_Force_vec;
output.SpikeTrain = spike_train;
output.Lce = Lce_vec;
output.Vce = Vce_vec;
output.Lse = Lse_vec;

%%
    function [twitch,T1,T2_temp] = twitch_function(Af,Lce,CT,RT,Fs)
        T1 = CT*Lce^2+CT*Af;
        T2_temp = (RT + RT*Af)/Lce;
        T2 = T2_temp/1.68;
        t_twitch = 0:1/Fs:2;
        f_1 = t_twitch./T1.*exp(1-t_twitch./T1);
        f_2 = t_twitch./T2.*exp(1-t_twitch./T2);
        
        twitch = [f_1(1:round(T1*Fs+1)) f_2(round(T2*Fs+1):end)];
        twitch = twitch(1:1*Fs);
        
    end

    function Af = Af_slow_function(f_env,L,Y)
        a_f = 0.56;
        n_f0 = 2.1;
        n_f1 = 5;
        n_f = n_f0 + n_f1*(1/L-1);
        Af = 1 - exp(-(Y*f_env/(a_f*n_f))^n_f);
    end

    function Af = Af_fast_function(f_env,L,S)
        a_f = 0.56;
        n_f0 = 2.1;
        n_f1 = 3.3;
        n_f = n_f0 + n_f1*(1/L-1);
        Af = 1 - exp(-(S*f_env/(a_f*n_f))^n_f);
        
    end

    function FF = Af_slow_correction_function(f_env,L,Y)
        a_f = 0.56;
        n_f0 = 2.1;
        n_f1 = 5;
        n_f = n_f0 + n_f1*(1/L-1);
        FF = 1 - exp(-(Y*f_env/(a_f*n_f))^n_f);
        offset = 0.375*L - 0.1775;
        alpha_FF = 1-offset;
        FF = FF*alpha_FF+offset;
        FF = FF/f_env;
        if FF > 1
            FF = 1;
        end
    end

    function FF = Af_fast_correction_function(f_env,L,S)
        a_f = 0.56;
        n_f0 = 2.1;
        n_f1 = 3.3;
        n_f = n_f0 + n_f1*(1/L-1);
        FF = 1 - exp(-(S*f_env/(a_f*n_f))^n_f);
        offset = 0.375*L - 0.1775;
        alpha_FF = 1-offset;
        FF = FF*alpha_FF+offset;
        FF = FF/f_env;
        if FF > 1
            FF = 1;
        end
    end

    function Y = yield_function(Y,Vce,Fs)
        c_y = 0.35;
        V_y = 0.1;
        T_y = 0.2;
        
        Y_dot = (1-c_y*(1-exp(-abs(Vce)/V_y))-Y)/T_y;
        Y = Y_dot*1/Fs + Y;
    end

    function S = sag_function(S,f_eff,Fs)
        if f_eff < 0.1
            a_s = 1.76;
        else
            a_s = 0.96;
        end
        T_s = 0.043;
        S_dot = (a_s-S)/T_s;
        S = S_dot*1/Fs + S;
    end

    function FL = FL_slow_function(L)
        %---------------------------
        % force length (F-L) relationship for slow-tiwtch fiber
        % input: normalized muscle length and velocity
        % output: F-L factor (0-1)
        %---------------------------
        beta = 2.3;
        omega = 1.12;
        rho = 1.62;
        
        FL = exp(-abs((L^beta - 1)/omega)^rho);
    end

    function FL = FL_fast_function(L)
        %---------------------------
        % force length (F-L) relationship for fast-twitch fiber
        % input: normalized muscle length and velocity
        % output: F-L factor (0-1)
        %---------------------------
        beta = 1.55;
        omega = 0.75;
        rho = 2.12;
        
        FL = exp(-abs((L^beta - 1)/omega)^rho);
    end

    function FVcon = FV_con_slow_function(L,V)
        %---------------------------
        % concentric force velocity (F-V) relationship for slow-twitch fiber
        % input: normalized muscle length and velocity
        % output: F-V factor (0-1)
        %---------------------------
        Vmax = -7.88;
        cv0 = 5.88;
        cv1 = 0;
        
        FVcon = (Vmax - V)/(Vmax + (cv0 + cv1*L)*V);
    end

    function FVcon = FV_con_fast_function(L,V)
        %---------------------------
        % concentric force velocity (F-V) relationship for fast-twitch fiber
        % input: normalized muscle length and velocity
        % output: F-V factor (0-1)
        %---------------------------
        Vmax = -9.15;
        cv0 = -5.7;
        cv1 = 9.18;
        
        FVcon = (Vmax - V)/(Vmax + (cv0 + cv1*L)*V);
    end

    function FVecc = FV_ecc_slow_function(L,V)
        %---------------------------
        % eccentric force velocity (F-V) relationship for slow-twitch fiber
        % input: normalized muscle length and velocity
        % output: F-V factor (0-1)
        %---------------------------
        av0 = -4.7;
        av1 = 8.41;
        av2 = -5.34;
        bv = 0.35;
        FVecc = (bv - (av0 + av1*L + av2*L^2)*V)/(bv+V);
    end

    function FVecc = FV_ecc_fast_function(L,V)
        %---------------------------
        % eccentric force velocity (F-V) relationship for fast-twitch fiber
        % input: normalized muscle length and velocity
        % output: F-V factor (0-1)
        %---------------------------
        av0 = -1.53;
        av1 = 0;
        av2 = 0;
        bv = 0.69;
        FVecc = (bv - (av0 + av1*L + av2*L^2)*V)/(bv+V);
    end

    function Fpe1 = F_pe_1_function(L,V)
        %---------------------------
        % passive element 1
        % input: normalized muscle length
        % output: passive element force (0-1)
        %---------------------------
        c1_pe1 = 23;
        k1_pe1 = 0.046;
        Lr1_pe1 = 1.17;
        eta = 0.01;
        
        Fpe1 = c1_pe1 * k1_pe1 * log(exp((L - Lr1_pe1)/k1_pe1)+1) + eta*V;
        
    end

    function Fpe2 = F_pe_2_function(L)
        %---------------------------
        % passive element 2
        % input: normalized muscle length
        % output: passive element force (0-1)
        %---------------------------
        c2_pe2 = -0.02;
        k2_pe2 = -21;
        Lr2_pe2 = 0.70;
        
        Fpe2 = c2_pe2*exp((k2_pe2*(L-Lr2_pe2))-1);
        
    end

    function Fse = F_se_function(LT)
        %---------------------------
        % series elastic element (tendon)
        % input: tendon length
        % output: tendon force (0-1)
        %---------------------------
        cT_se = 27.8; %27.8
        kT_se = 0.0047;
        LrT_se = 0.964;
        
        Fse = cT_se * kT_se * log(exp((LT - LrT_se)/kT_se)+1);
        
    end

    function [Lce_initial,Lse_initial,Lmax] =  InitialLength_function(L0,alpha,Lse_slack,Lce_initial,Lse_initial)
        %---------------------------
        % Determine the initial lengths of muscle and tendon and maximal
        % muscle length
        %---------------------------
        
        % serires elastic element parameters
        cT = 27.8;
        kT = 0.0047;
        LrT = 0.964;
        % parallel passive element parameters
        c1 = 23;
        k1 = 0.046;
        Lr1 = 1.17;
        
        % passive force produced by parallel passive element at maximal
        % muscle length
        PassiveForce = c1 * k1 * log(exp((1 - Lr1)/k1)+1);
        % tendon length at the above passive force
        Normalized_SE_Length = kT*log(exp(PassiveForce/cT/kT)-1)+LrT;
        
        % maximal musculotendon length defined by joint range of motion
        Lmt_temp_max = L0*cos(alpha)+Lse_slack + 1;
        
        % optimal muscle length
        L0_temp = L0;
        % optimal tendon length (Song et al. 2008)
        L0T_temp = Lse_slack*1.05;
        
        % tendon length at maximal muscle length
        SE_Length =  L0T_temp * Normalized_SE_Length;
        % maximal fasicle length
        FasclMax = (Lmt_temp_max - SE_Length)/L0_temp;
        % maximal muscle fiber length
        Lmax = FasclMax/cos(alpha);
        
        % initial musculotendon length defined by the user input
        Lmt_temp = Lce_initial * cos(alpha) + Lse_initial;
        
        % initial muscle length determined by passive muscle force and
        % tendon force
        InitialLength =  (Lmt_temp-(-L0T_temp*(kT/k1*Lr1-LrT-kT*log(c1/cT*k1/kT))))/(100*(1+kT/k1*L0T_temp/Lmax*1/L0_temp)*cos(alpha));
        % normalize the muscle legnth to optimal muscle length
        Lce_initial = InitialLength/(L0_temp/100);
        % calculate initial length of tendon and normalize it to optimal
        % tendon length
        Lse_initial = (Lmt_temp - InitialLength*cos(alpha)*100)/L0T_temp;
    end

end