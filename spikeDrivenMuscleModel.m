%==========================================================================
% spikeDrivenMuscleModel.m
% Author: Akira Nagamori
% Last update: 3/5/19
%==========================================================================
function [output] = spikeDrivenMuscleModel(Fs,time,input,modelParameter,figOpt)
%% Simulation parameters
synaptic_drive = input;

%%
recruitmentType = modelParameter.recruitment;

%% Muscle architectural parameters
L0 = modelParameter.optimalLength; % optimal muscle length [cm]
mass = modelParameter.mass; % muscle mass [kg]
F0 = modelParameter.F0; % maximal force

L0T = modelParameter.L0T;
alpha = modelParameter.pennationAngle;
Lmt =modelParameter.Lmt; % intial musculotendon length
L_ce = modelParameter.L_ce;
L_se = modelParameter.L_se;
Lmax = modelParameter.Lmax;
%% Motor unit architecture
N_MU = modelParameter.N_MU; % number of motor units
i_MU = modelParameter.i_MU; % index for motor units

%% Peak tetanic force
PTi_new = modelParameter.PTi_new;

%% Fractional PSCA
%F_pcsa_slow = 0.3; % fractional PSCA of slow-twitch motor units (0-1)
%[~, index_slow] = min(abs(cumsum(PTi) - F0*F_pcsa_slow)); rng(1)
index_slow = modelParameter.index_slow;

%% Motor unit parameter
parameter_Matrix = modelParameter.parameterMatrix;

%% Recruitment threshold
U_th_new = modelParameter.U_th_new;

%% Minimum and maximum firing rate
FR_half = modelParameter.FR_half;
MDR = modelParameter.MDR;
PDR = modelParameter.PDR;

g_e = modelParameter.g_e;
if recruitmentType == 3
    index_saturation = modelParameter.index_saturation;
    lamda = modelParameter.lamda;
    k_e = modelParameter.k_e;
    U_th_t = modelParameter.U_th_t;
end
%% Activation dynamics (Song et al., 2008)
U_eff = 0;
T_U = 0.03;

tau_1 = parameter_Matrix(:,9);
R_temp = exp(-time./tau_1);
gamma = parameter_Matrix(:,15);

%% Sag parameter
a_s = ones(N_MU,1)*0.96;
%% Discharge rate parameters
cv_MU = modelParameter.CV_MU; % coefficient of variation for interspike intervals

%% Muscle length
V_ce = 0;
%% Initilization
DR_temp = zeros(N_MU,1);
DR_mat = zeros(N_MU,length(time));

spike_time = zeros(N_MU,1);
spike_train = zeros(N_MU,length(time));
force = zeros(N_MU,length(time));
F_ce = zeros(1,length(time));
F_se = zeros(1,length(time));
F_total = zeros(1,length(time));

R = zeros(N_MU,length(time));
c = zeros(N_MU,1);
cf = zeros(N_MU,1);
A = zeros(N_MU,1);
c_mat = zeros(N_MU,length(time));
cf_mat = zeros(N_MU,length(time));
A_tilde_mat = zeros(N_MU,length(time));
A_mat = zeros(N_MU,length(time));

S_i = zeros(N_MU,1);
Y_i = zeros(N_MU,1);
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
            DR_temp = g_e.*(U_eff-U_th_new)+0.5;
            DR_temp(DR_temp<0.5) = 0;
            DR_temp(DR_temp>2) = 2;
            DR_MU = DR_temp.*FR_half;
        elseif recruitmentType == 2
            DR_MU = g_e.*(U_eff-U_th_new)+MDR;
            DR_MU(DR_MU<MDR) = 0;
            DR_MU(DR_MU>PDR) = PDR(DR_MU>PDR);
        elseif recruitmentType == 3
            DR_MU = g_e.*(U_eff-U_th_new)+MDR;
            for m = 1:index_slow
                index = index_saturation(m);
                if U_eff <= U_th_t(index)
                    DR_temp(index) = MDR(index) + lamda(index).*k_e(index)*(U_eff-U_th_new(index));
                else
                    DR_temp(index) = PDR(index)-k_e(index)*(1-U_eff);
                end
            end
            DR_MU(index_saturation) = DR_temp(index_saturation);
            DR_MU(DR_MU<MDR) = 0;
            DR_MU(DR_MU>PDR) = PDR(DR_MU>PDR);
        end
        % Zero the discharge rate of a MU if it is smaller than its minimum
        % firing rate
        
        DR_mat(:,t) = DR_MU;
        %% Sag & Yield (Song et al., 2008)
        f_eff = DR_MU./FR_half;
        S_i = sag_function(S_i,f_eff,a_s,Fs);
        S_i(1:index_slow) = 1;
        S_mat(:,t) = S_i;
        Y_i = yield_function(Y_i,V_ce,Fs);
        Y_i(index_slow+1:end) = 1;
        Y_mat(:,t) = Y_i;
        
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
                
                temp = conv(spike_train_temp,R_temp(n,:)*(1+2*A(n)^gamma(n)));
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
                    
                    temp = conv(spike_train_temp,R_temp(n,:)*(1+2*A(n)^gamma(n)));
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
                    
                    temp = conv(spike_train_temp,R_temp(n,:)*(1+2*A(n)^gamma(n)));
                    R(n,:) = R(n,:) + temp(1:length(time));
                end
            end
        end
    end
    %% Convert spikes into activation
    [c,cf,A_tilde,A] = spike2activation(R(:,t),c,cf,A,parameter_Matrix,L_ce,S_i,Y_i,Fs);
    
    c_mat(:,t) = c;
    cf_mat(:,t) = cf;
    A_tilde_mat(:,t) = A_tilde;
    A_mat(:,t) = A;
    
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
    
    f_i = A.*PTi_new'.*(FL.*FV+F_pe2);
    %f_i = A.*PTi_new'.*(FL+F_pe2);
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
if figOpt == 1
    figure(1)
    plot(time,F_ce)
    xlabel('Time (s)')
    ylabel('Force (N)')
    hold on
end


output.spike_train = spike_train;
output.Force = F_ce;
output.force = force;
output.ForceTendon = F_se;
output.Lce = MuscleLength./(L0/100);
output.Vce = MuscleVelocity./(L0/100);

%% Convert spike trian into activation
    function [c,cf,A_tilde,A] = spike2activation(R,c,cf,A,parameter_Matrix,Lce,S_i,Y_i,Fs)
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
        
        %%
        c_dot = k_1.*(C-c-cf).*R - k_2.*c.*(S-C+c+cf)-(k_3.*c-k_4.*cf).*(1-cf);
        cf_dot = (1-cf).*(k_3.*c-k_4.*cf);
        c = c_dot/Fs + c;
        cf = cf_dot/Fs + cf;
        
        %% Stage 2
        % Cooperativity and saturation
        if cf < 0
            cf_temp = 0;
        else
            cf_temp = cf.*S_i.*Y_i;
        end
        A_tilde = cf_temp.^N./(cf_temp.^N+K.^N);
        
        %% Stage 3
        % First-order dynamics to muscle activation, A
        %tau_2 = tau_2_0*(1-0.8*A_tilda)^2;
        A_dot = (A_tilde-A)./tau_2;
        A = A_dot./Fs + A;
        
    end

%% Sag
    function [S] = sag_function(S,f_eff,a_s,Fs)
        
        a_s(f_eff<0.1) = 20;
        
        T_s = 0.015;
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

end