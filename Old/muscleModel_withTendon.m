%==========================================================================
% muscleModel_withTendon.m
% Author: Akira Nagamori
% Last update: 3/5/19
% Descriptions:
%   Full model without tendon
%==========================================================================
function [output] = muscleModel_withTendon(Fs,time,input,modelParameter)
%% Simulation parameters
synaptic_drive = input;

%% 
recruitmentType = modelParameter.recruitment;

%% Muscle architectural parameters
density = 1.06; %
L0 = modelParameter.optimalLength; % optimal muscle length [cm]
mass = modelParameter.mass; % muscle mass [kg]
PCSA = (mass*1000)/(density*L0); % PCSA of muscle
sigma = 31.8; % specific tension
F0 = PCSA * sigma; % maximal force

L_tendon = modelParameter.tendonSlackLength;
L0T = L_tendon*1.05;
alpha = modelParameter.pennationAngle;
Lm_initial = modelParameter.muscleInitialLength; % muscle initial length
Lt_initial = modelParameter.tendonInitialLength; % tendon initial length
Lmt = Lm_initial*cos(alpha)+Lt_initial; % intial musculotendon length
[L_ce,L_se,Lmax] =  InitialLength_function(modelParameter);

%--------------------------------------------------------------------------
% Motor unit architecture
N_MU = 300; % number of motor units
i_MU = 1:N_MU; % index for motor units

%% Peak tetanic force
RP_MU = modelParameter.RP; %range of peak tension across motor untis in unit of fold
b_MU = log(RP_MU)/N_MU; %coefficient to establish a range of twich force values
P_MU = exp(b_MU*i_MU); %force generated by a motor unit as a function of its recruitment threshold
PTi = P_MU./sum(P_MU)*F0; % peak tetanic force for individual units

%% Fractional PSCA
F_pcsa_slow = 0.3; % fractional PSCA of slow-twitch motor units (0-1)
[~, index_slow] = min(abs(cumsum(PTi) - F0*F_pcsa_slow)); rng(1)
%index_slow = 196;
%% Motor unit parameter
parameter_Matrix = modelParameter.parameterMatrix;
%matrix of [N_MU,15]
% N_MU is in order of CT_sorted

%% Assign peak tetanic force into each unit
% shuffle indexes within each fiber type with respect to contraction time
% this will allow us to randomly assign peak tetanic tension to each motor
% unit with different contraction time
rng(1)
R_slow = randperm(index_slow);
index_fast = index_slow+1:N_MU;
R_fast_temp = randperm(length(index_fast));
R_fast = index_fast(R_fast_temp);
index_MU_PTi = [R_slow R_fast]; % vector of indexes to match peak tetanic tension to appropriate contraction time
PTi_new = PTi (index_MU_PTi);

%% Recruitment threshold
% Find recruitment threshold for individual units using exponential fit
% Recruitment threshold is correlated to peak tetanic tension
%   Use index_MU_PTi to appropriately index each MU
Ur = modelParameter.Ur; % recruitment threshold for the lastly recruited motor unit
Ur_1 = 0.01; % reruitment threshold for the first unit
f_RT = fit([1 N_MU]',[Ur_1 Ur]','exp1');
coeffs_f_RT = coeffvalues(f_RT);
U_th = coeffs_f_RT(1)*exp(coeffs_f_RT(2)*i_MU); % the resulting recruitment threshold for individual units
U_th_new = U_th(index_MU_PTi);


%% Minimum and maximum firing rate
FR_half = modelParameter.FR_half;
MDR = FR_half/2;
PDR = FR_half*2;

%% Gain for frequency-activation relationship
g_e = max((PDR-MDR)./(1-U_th_new));

%% Activation dynamics (Song et al., 2008)
U_eff = 0;
T_U = 0.03;

tau_1 = parameter_Matrix(:,9);
R_temp = exp(-time./tau_1);
alpha_MU = parameter_Matrix(:,15);

%% Sag parameter
a_s = ones(N_MU,1)*0.96;
%% Discharge rate parameters
cv_MU = modelParameter.CV_MU; % coefficient of variation for interspike intervals

%% Muscle length
V_ce = 0;
%% Initilization
DR_mat = zeros(N_MU,length(time));

spike_time = zeros(N_MU,1);
spike_train = zeros(N_MU,length(time));
force = zeros(N_MU,length(time));
F_ce = zeros(1,length(time));
F_se = zeros(1,length(time));
F_total = zeros(1,length(time));

R = zeros(N_MU,length(time));
x = zeros(N_MU,1);
y = zeros(N_MU,1);
z = zeros(N_MU,1);

x_mat = zeros(N_MU,length(time));
y_mat = zeros(N_MU,length(time));
z_mat = zeros(N_MU,length(time));

S_MU = zeros(N_MU,1);
Y = zeros(N_MU,1);
S_mat = zeros(N_MU,length(time));
Y_mat = zeros(N_MU,length(time));

FL = zeros(N_MU,1);
FV = zeros(N_MU,1);

MuscleVelocity = zeros(1,length(time));
MuscleLength = zeros(1,length(time));
MuscleLength(1) = L_ce*L0/100;

h = 1/Fs;
%% Simulation
rng('shuffle')
for t = 1:length(time)
    if t > 1
        %% Effective activation (Song et al., 2008)
        U_eff_dot = (synaptic_drive(t) - U_eff)/T_U;
        U_eff = U_eff_dot*1/Fs + U_eff;
        
        %% Calculate firing rate
        % Linear increase in discharge rate up to Ur
        if recruitmentType == 1
            DR_MU = g_e.*(U_eff-U_th_new)+MDR;
        elseif recruitmentType == 2
            DR_MU = g_e.*(U_eff-U_th_new)+MDR;
        elseif recruitmentType == 3
            DR_MU= PDR.*(1-exp(-(U_eff-U_th_new)./g_e))+MDR;
        end
        % Zero the discharge rate of a MU if it is smaller than its minimum
        % firing rate
        DR_MU(DR_MU<MDR) = 0;
        DR_MU(DR_MU>PDR) = PDR(DR_MU>PDR);
        DR_mat(:,t) = DR_MU;
        %% Sag & Yield (Song et al., 2008)
        f_eff = DR_MU./FR_half;
        S_MU = sag_function(S_MU,f_eff,a_s,Fs);
        S_MU(1:index_slow) = 1;
        S_mat(:,t) = S_MU;
        Y = yield_function(Y,V_ce,Fs);
        Y(index_slow+1:end) = 1;
        Y_mat(:,t) = Y;
        
        %% Convert activation into spike trains
        index_1 = i_MU(DR_MU >= MDR & DR_mat(:,t-1)' == 0);
        index_2 = i_MU(DR_MU >= MDR & spike_time'==t);
        index = [index_1 index_2];
        
        for j = 1:length(index) % loop through motor units whose firing rate is greater than minimum firing rate defined by the user
            n = index(j);
            spike_train_temp = zeros(1,length(t));
            if ~any(spike_train(n,:)) % when the motor unit fires at the first time
                spike_train(n,t) = 1; % add a spike to the vector
                spike_train_temp(t) = 1;
                mu = 1/DR_MU(n);
                Z = randn(1);
                Z(Z>3.9) = 3.9;
                Z(Z<-3.9) = -3.9;                
                spike_time_temp = (mu + mu*cv_MU*Z)*Fs;
                if spike_time_temp <= 0.002*Fs
                    spike_time_temp = 0.002*Fs;
                end
                spike_time(n) = round(spike_time_temp) + t;
                
                temp = conv(spike_train_temp,R_temp(n,:)*(1+2*z(n)^alpha_MU(n)));
                R(n,:) = R(n,:) + temp(1:length(time));
            else % when the motor unit have already fired at least once
                if spike_time(n) == t % when the motor unit fires
                    spike_train(n,t) = 1;
                    spike_train_temp(t) = 1;
                    % update mean firing rate of the motor unit given the
                    % current value of input
                    mu = 1/DR_MU(n); % interspike interval
                    Z = randn(1);
                    Z(Z>3.9) = 3.9;
                    Z(Z<-3.9) = -3.9;
                    spike_time_temp = (mu + mu*cv_MU*Z)*Fs; % interspike interval
                    if spike_time_temp <= 0.002*Fs
                        spike_time_temp = 0.002*Fs;
                    end
                    spike_time(n) = round(spike_time_temp) + t;
                    
                    temp = conv(spike_train_temp,R_temp(n,:)*(1+2*z(n)^alpha_MU(n)));
                    R(n,:) = R(n,:) + temp(1:length(time));
                elseif t > spike_time(n) + round(1/DR_MU(n)*Fs)
                    spike_train(n,t) = 1;
                    spike_train_temp(t) = 1;
                    spike_time(n) = t;
                    mu = 1/DR_MU(n); % interspike interval
                    Z = randn(1);
                    Z(Z>3.9) = 3.9;
                    Z(Z<-3.9) = -3.9;
                    spike_time_temp = (mu + mu*cv_MU*Z)*Fs; % interspike interval
                    if spike_time_temp <= 0.002*Fs
                        spike_time_temp = 0.002*Fs;
                    end
                    spike_time(n) = round(spike_time_temp) + t;
                    
                    temp = conv(spike_train_temp,R_temp(n,:)*(1+2*z(n)^alpha_MU(n)));
                    R(n,:) = R(n,:) + temp(1:length(time));
                end
            end
        end
    end
        %% Convert spikes into activation
        [x,y,z] = spike2activation(R(:,t),x,y,z,parameter_Matrix,L_ce,Fs);
        
        x_mat(:,t) = x;
        y_mat(:,t) = y;
        z_mat(:,t) = z;
        
        %% Force-length and force-velocity
        FL(1:index_slow) = FL_slow_function(L_ce);
        FL(index_slow+1:end) = FL_fast_function(L_ce);
        
        if V_ce > 0
            FV(1:index_slow) = FVecc_slow_function(L_ce,V_ce);
            FV(index_slow+1:end) = FVecc_fast_function(L_ce,V_ce);
        else
            FV(1:index_slow) = FVcon_slow_function(L_ce,V_ce);
            FV(index_slow+1:end) = FVcon_fast_function(L_ce,V_ce);
        end
        
        %% Passive element 1
        F_pe1 = Fpe1_function(L_ce/Lmax,V_ce);
        
        %% Passive element 2
        F_pe2 = Fpe2_function(L_ce);
        if F_pe2 > 0
            F_pe2 = 0;
        end
        
        %f_i = z.*PTi_new'.*S_MU.*Y.*(FL.*FV+F_pe2);
        f_i = z.*PTi_new'.*S_MU.*Y.*(FL+F_pe2);
        force(:,t) = f_i;
        
        F_ce(t) = sum(f_i);
        F_total(t) = F_ce(t) + F_pe1*F0;
        
        F_se(t) = Fse_function(L_se) * F0;
        
        k_0_de = h*MuscleVelocity(t);
        l_0_de = h*((F_se(t)*cos(alpha) - F_total(t)*(cos(alpha)).^2)/(mass) ...
            + (MuscleVelocity(t)).^2*tan(alpha).^2/(MuscleLength(t)));
        k_1_de = h*(MuscleVelocity(t)+l_0_de/2);
        l_1_de = h*((F_se(t)*cos(alpha) - F_total(t)*(cos(alpha)).^2)/(mass) ...
            + (MuscleVelocity(t)+l_0_de/2).^2*tan(alpha).^2/(MuscleLength(t)+k_0_de/2));
        k_2_de = h*(MuscleVelocity(t)+l_1_de/2);
        l_2_de = h*((F_se(t)*cos(alpha) - F_total(t)*(cos(alpha)).^2)/(mass) ...
            + (MuscleVelocity(t)+l_1_de/2).^2*tan(alpha).^2/(MuscleLength(t)+k_1_de/2));
        k_3_de = h*(MuscleVelocity(t)+l_2_de);
        l_3_de = h*((F_se(t)*cos(alpha) - F_total(t)*(cos(alpha)).^2)/(mass) ...
            + (MuscleVelocity(t)+l_2_de).^2*tan(alpha).^2/(MuscleLength(t)+k_2_de));
        MuscleLength(t+1) = MuscleLength(t) + 1/6*(k_0_de+2*k_1_de+2*k_2_de+k_3_de);
        MuscleVelocity(t+1) = MuscleVelocity(t) + 1/6*(l_0_de+2*l_1_de+2*l_2_de+l_3_de);
        
        % normalize each variable to optimal muscle length or tendon length
        V_ce = MuscleVelocity(t+1)/(L0/100);
        L_ce = MuscleLength(t+1)/(L0/100);
        L_se = (Lmt - L_ce*L0*cos(alpha))/L0T;
end

%%
figure(1)
plot(time,F_ce)
xlabel('Time (s)')
ylabel('Force (N)')
hold on


output.spike_train = spike_train;
output.Force = F_ce;
output.force = force;
output.ForceTendon = F_se;
output.Lce = MuscleLength./(L0/100);
output.Vce = MuscleVelocity./(L0/100);

%% Convert spike trian into activation
    function [x,y,z] = spike2activation(R,x,y,z,parameter_Matrix,Lce,Fs)
        S = parameter_Matrix(:,1); %7;
        C = parameter_Matrix(:,2); %1.025;
        k_1 = parameter_Matrix(:,3); %14.625;
        k_2 = parameter_Matrix(:,4); %4.9375;
        k_3 = parameter_Matrix(:,5)*Lce + parameter_Matrix(:,6); %17.41*Lce - 2.85;
        k_4 = parameter_Matrix(:,7)*Lce + parameter_Matrix(:,8); %-7.67*Lce + 14.92;
        tau_2 = parameter_Matrix(:,10); % 0.04;
        N = parameter_Matrix(:,11)*Lce + parameter_Matrix(:,12); %-2.26*Lce + 4.20;
        K = parameter_Matrix(:,13)*Lce + parameter_Matrix(:,14); %-0.044*Lce + 0.080;
        %4.475;
        
        x_dot = k_1.*(C-x-y).*R - k_2.*x.*((C.*(S-1))+x+y)-(k_3.*x+k_4.*y).*(1-y);
        y_dot = (1-y).*(k_3.*x-k_4.*y);
        x = x_dot/Fs + x;
        y = y_dot/Fs + y;
        y_temp = y;
        y_temp(y_temp<0) = 0;
        
        y_int = y_temp.^N./(y_temp.^N+K.^N);
        
        z_dot = (y_int-z)./tau_2;
        z = z_dot/Fs + z; % activation
        
    end

%% Sag
    function [S] = sag_function(S,f_eff,a_s,Fs)
        
        a_s(f_eff<0.1) = 1.76;
        
        T_s = 0.043;
        S_dot = (a_s - S)./T_s;
        S = S_dot/Fs + S;
        
    end

%% Yield
    function [Y] = yield_function(Y,V,Fs)
        c_y = 0.35;
        V_y = 0.1;
        T_y = 0.2;
        Y_dot = (1-c_y.*(1-exp(-abs(V)./V_y))-Y)./T_y;
        Y = Y_dot/Fs + Y;
        
    end

%% Force-length relationship for slow twitch
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

%% Force-length relationship for fast twitch
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

%% Concentric force-velocity relationship for slow twitch
    function FVcon = FVcon_slow_function(L,V)
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

%% Concentric force-velocity relationship for fast twitch
    function FVcon = FVcon_fast_function(L,V)
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

%% Eccentric force-velocity relationship for slow twitch
    function FVecc = FVecc_slow_function(L,V)
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

%% Eccentric force-velocity relationship for slow twitch
    function FVecc = FVecc_fast_function(L,V)
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

%% Force-length relationship for passive element 1
    function Fpe1 = Fpe1_function(L,V)
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

%% Force-length relationship for passive element 2
    function Fpe2 = Fpe2_function(L)
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

%% Force-length relationship for series-elastic element
    function Fse = Fse_function(LT)
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

    function [Lce_initial,Lse_initial,Lmax] =  InitialLength_function(modeParameter)
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
        Lmt_temp_max = modeParameter.optimalLength*cos(modeParameter.pennationAngle) ...
            +modeParameter.tendonSlackLength + 1;
        
        % optimal muscle length
        L0_temp = modeParameter.optimalLength;
        % optimal tendon length (Song et al. 2008)
        L0T_temp = modeParameter.tendonSlackLength*1.05;
        
        % tendon length at maximal muscle length
        SE_Length =  L0T_temp * Normalized_SE_Length;
        % maximal fasicle length
        FasclMax = (Lmt_temp_max - SE_Length)/L0_temp;
        % maximal muscle fiber length
        Lmax = FasclMax/cos(modeParameter.pennationAngle);
        
        % initial musculotendon length defined by the user input
        Lmt_temp = modeParameter.muscleInitialLength * cos(modeParameter.pennationAngle) + modeParameter.tendonInitialLength;
        
        % initial muscle length determined by passive muscle force and
        % tendon force
        InitialLength =  (Lmt_temp-(-L0T_temp*(kT/k1*Lr1-LrT-kT*log(c1/cT*k1/kT))))/(100*(1+kT/k1*L0T_temp/Lmax*1/L0_temp)*cos(modeParameter.pennationAngle));
        % normalize the muscle legnth to optimal muscle length
        Lce_initial = InitialLength/(L0_temp/100);
        % calculate initial length of tendon and normalize it to optimal
        % tendon length
        Lse_initial = (Lmt_temp - InitialLength*cos(modeParameter.pennationAngle)*100)/L0T_temp;
    end

end