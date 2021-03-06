%==========================================================================
% spikeDrivenMuscleModel_MN_SLR.m
% Author: Akira Nagamori
% Last update: 8/19/2019
%==========================================================================
function [output] = spikeDrivenMuscleModel_MN_SLR(Fs,time,input,modelParameter,parameterMN,SLRParameter,controlOpt,figOpt)
%% Simulation parameters
C_input = input;

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

%% Peak tetanic force
PTi = modelParameter.PTi;

%% Fractional PSCA
%F_pcsa_slow = 0.3; % fractional PSCA of slow-twitch motor units (0-1)
%[~, index_slow] = min(abs(cumsum(PTi) - F0*F_pcsa_slow)); rng(1)
index_slow = modelParameter.index_slow;

%% Motor unit parameter
parameter_Matrix = modelParameter.parameterMatrix;

%% Recruitment threshold
U_th = modelParameter.U_th;

%% Minimum and maximum firing rate
FR_half = modelParameter.FR_half;

g_e = modelParameter.g_e;

if recruitmentType == 3
    index_t  = modelParameter.index_t; % row vector
    lamda = modelParameter.lamda; % row vector
    k_e = modelParameter.k_e; % row vector
    U_th_t = modelParameter.U_th_t; % row vector
end
%% Activation dynamics (Song et al., 2008)
tau_1 = parameter_Matrix(:,9);
R_temp = exp(-time./tau_1);
gamma = parameter_Matrix(:,15);

%% Sag parameter
a_s = ones(N_MU,1)*0.96;

%% Motoneuron parameters
a_MN = parameterMN.a;
I_th = parameterMN.I_th;
I_max = parameterMN.I_max;

%% Muscle length
A_ce = 0;
V_ce = 0;
%% Initilization

force = zeros(N_MU,length(time));
F_ce = zeros(1,length(time));
F_se = zeros(1,length(time));
F_total = zeros(1,length(time));

% motoneuron related parameters
DR_MU = zeros(N_MU,1);
spike_train = zeros(N_MU,length(time));

I_initial = (I_max-I_th)./(1-U_th).*(0-U_th) + I_th;
v_MN = (-5+0.2-sqrt((5-0.2)^2-4*0.04*(140+I_initial)))/(2*0.04);
u_MN = 0.2*v_MN;
%v_MN = ones(N_MU,1)*-65;
v_mat = zeros(N_MU,length(time));
I_mat = zeros(N_MU,length(time));

% activation-frequency
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
MuscleAcceleration = zeros(1,length(time));
MuscleLength(1) = L_ce*L0/100;

%%
gamma_dynamic = SLRParameter.gamma_dynamic;
gamma_static = SLRParameter.gamma_static;
Ia_delay = SLRParameter.Ia_delay;
Ia_gain = SLRParameter.Ia_gain;

f_dynamic_bag1 = 0;
T_bag1 = 0;
T_dot_bag1 = 0;
f_static_bag2 = 0;
T_bag2 = 0;
T_dot_bag2 = 0;
T_chain = 0;
T_dot_chain = 0;

FR_Ia = zeros(1,length(time));
Ia_Input = zeros(N_MU,length(time));

%%
Ib_gain = SLRParameter.Ib_gain;
Ib_delay = SLRParameter.Ib_delay;
G1 = SLRParameter.G1;
G2 = SLRParameter.G2;
num1_GTO = SLRParameter.num_GTO(1);
num2_GTO = SLRParameter.num_GTO(2);
num3_GTO = SLRParameter.num_GTO(3);
den1_GTO = SLRParameter.den_GTO(1);
den2_GTO = SLRParameter.den_GTO(2);
den3_GTO = SLRParameter.den_GTO(3);

x_GTO = zeros(1,length(time));
FR_Ib_temp = zeros(1,length(time));
FR_Ib = zeros(1,length(time));
Ib_Input = zeros(N_MU,length(time));

%%
RI_gain = SLRParameter.Ib_gain;
RI_delay = SLRParameter.RI_delay;
num1_RI = SLRParameter.num_RI(1);
num2_RI = SLRParameter.num_RI(2);
num3_RI = SLRParameter.num_RI(3);
den1_RI = SLRParameter.den_RI(1);
den2_RI = SLRParameter.den_RI(2);
den3_RI = SLRParameter.den_RI(3);

FR_RI_temp = zeros(N_MU,length(time));
FR_RI = zeros(N_MU,length(time));
RI_Input = zeros(N_MU,length(time));

%%
C_delay = SLRParameter.C_delay;
K_C = SLRParameter.K_C;
%%
noise_amp_Ia = SLRParameter.noise_amp_Ia;
noise_amp_Ib = SLRParameter.noise_amp_Ib;
noise_amp_RI = SLRParameter.noise_amp_RI;
noise_amp_C = SLRParameter.noise_amp_C;
noise_amp_ID = SLRParameter.noise_amp_ID;
noise_amp_CD = SLRParameter.noise_amp_CD;
noise_Ia =  zeros(N_MU,1);
noise_Ib =  zeros(N_MU,1);
noise_RI =  zeros(N_MU,1);
noise_C =  zeros(N_MU,1);
noise_ID = zeros(N_MU,1);
noise_CM = zeros(1,1);
%%
C_temp = 0;
U_mat = zeros(N_MU,length(time));
%%
h = 1/Fs;
%% Simulation
rng('shuffle')
for t = 1:length(time)
    %%
    [AP_bag1,f_dynamic_bag1,T_bag1,T_dot_bag1] = bag1_model(f_dynamic_bag1,gamma_dynamic,T_bag1,T_dot_bag1,L_ce,V_ce,A_ce,Fs);
    [AP_primary_bag2,AP_secondary_bag2,f_static_bag2,T_bag2,T_dot_bag2] = bag2_model(f_static_bag2,gamma_static,T_bag2,T_dot_bag2,L_ce,V_ce,A_ce,Fs);
    [AP_primary_chain,AP_secondary_chain,T_chain,T_dot_chain] = chain_model(gamma_static,T_chain,T_dot_chain,L_ce,V_ce,A_ce,Fs);
    [Output_Primary,~] = SpindleOutput(AP_bag1,AP_primary_bag2,AP_secondary_bag2,AP_primary_chain,AP_secondary_chain);
    
    FR_Ia(t) = Output_Primary;
    Ia_Input(:,t) = FR_Ia(t)/Ia_gain;
    [noise_Ia] = noise(noise_Ia,noise_amp_Ia,Fs);
    Ia_Input(:,t) = Ia_Input(:,t) + Ia_Input(:,t).*noise_Ia;
    %%
    if t > 5
        %% GTO activity
        x_GTO(t) = G1*log(F_se(t-1)/G2+1);
        FR_Ib_temp(t) = (num3_GTO*x_GTO(t-2) + num2_GTO*x_GTO(t-1) + num1_GTO*x_GTO(t) - ...
            den3_GTO*FR_Ib_temp(t-2) - den2_GTO*FR_Ib_temp(t-1))/den1_GTO;
        FR_Ib(t) = FR_Ib_temp(t);
        if FR_Ib(t)<0
            FR_Ib(t) = 0;
        end
        %% Renshaw cell activity
        FR_RI_temp(:,t) = (num3_RI*U_mat(:,t-2) + num2_RI*U_mat(:,t-1) + num1_RI*U_mat(:,t)...
            - den3_RI*FR_RI_temp(:,t-2) - den2_RI*FR_RI_temp(:,t-1))/den1_RI;
        FR_RI(:,t) = FR_RI_temp(:,t);
        index_negative = FR_RI(:,t) < 0;
        FR_RI(index_negative,t) = 0;
    end
    Ib_Input(:,t) = FR_Ib(t)/Ib_gain;
    [noise_Ib] = noise(noise_Ib,noise_amp_Ib,Fs);
    Ib_Input(:,t) = Ib_Input(:,t) + Ib_Input(:,t).*noise_Ib;
    
    RI_Input(:,t) = FR_RI(:,t)/RI_gain;
    [noise_RI] = noise(noise_RI,noise_amp_RI,Fs);
    RI_Input(:,t) = RI_Input(:,t) + RI_Input(:,t).*noise_RI;
    
    %%
    [noise_C] = noise(noise_C,noise_amp_C,Fs);
    [noise_ID] = noise(noise_ID,noise_amp_ID,Fs);
    [noise_CM] = noise(noise_CM,noise_amp_CD,Fs);
    %%
    if t > 1
        %% Control input, U, to the entire pool
        if controlOpt == 1
            % Feedfoward input
            U = C_input(t);
        elseif controlOpt == 2
            % Feedback input
            if t > Ia_delay && t <= Ib_delay
                U = C_input(t) + noise_C*C_input(t) ...
                    + Ia_Input(:,t-Ia_delay) ...
                    - RI_Input(:,t-RI_delay);
            elseif t > Ib_delay && t <= C_delay
                U = C_input(t) + noise_C*C_input(t) ...
                    + Ia_Input(:,t-Ia_delay) ...
                    - Ib_Input(:,t-Ib_delay)...
                    - RI_Input(:,t-RI_delay);
            elseif t > C_delay
                C_temp = K_C*(C_input(t) - F_se(t-C_delay)/F0) + C_temp;
                U = C_temp + noise_C*C_temp...
                    + Ia_Input(:,t-Ia_delay) ...
                    - Ib_Input(:,t-Ib_delay)...
                    - RI_Input(:,t-RI_delay);
            else
                U = C_input(t) + noise_C*C_input(t);
            end
            if U < 0
                U = 0;
            end
        end
        U = U; %+noise_CM*U;
        U_mat(:,t) = U;
        
        %% Calculate firing rate
        % Linear increase in discharge rate up to Ur
        if recruitmentType == 1
            I = g_e.*(U+noise_ID.*U-U_th) + I_th;
        elseif recruitmentType == 3
            I = zeros(N_MU,1);
            U_temp = U; % + noise_ID*U;
            I_temp_1 = I_th + lamda.*k_e.*(U_temp-U_th);
            index_1 = find(U_temp <= U_th_t);
            I(index_1) = I_temp_1(index_1);
            I_temp_2 = I_max-k_e.*(1-U_temp);
            index_2 = find(U_temp > U_th_t);
            I(index_2) = I_temp_2(index_2);
            I_temp_3 = g_e.*(U_temp-U_th)+I_th;
            I(index_t) = I_temp_3(index_t);
            I = (I+I.*noise_ID+I*noise_CM);
        end
        % Zero the discharge rate of a MU if it is smaller than its minimum
        % firing rate
        spike_vec = zeros(N_MU,1);
        [u_MN,v_MN,spike_vec] = motoneuron(I,u_MN,v_MN,spike_vec,a_MN,Fs);
        spike_train(:,t) = spike_vec;
        v_mat(:,t) = v_MN;
        I_mat(:,t) = I;
        
        %% Sag & Yield (Song et al., 2008)
        f_eff = DR_MU./FR_half;
        S_i = sag_function(S_i,f_eff,a_s,Fs);
        S_i(1:index_slow) = 1;
        S_mat(:,t) = S_i;
        Y_i = yield_function(Y_i,V_ce,Fs);
        Y_i(index_slow+1:end) = 1;
        Y_mat(:,t) = Y_i;
        
        %% Convert activation into spike trains
        index =  find(spike_vec==1);
        
        for j = 1:length(index) % loop through motor units whose firing rate is greater than minimum firing rate defined by the user
            n = index(j);
            spike_train_temp = zeros(1,length(t));
            spike_train_temp(t) = 1;
            
            temp = conv(spike_train_temp,R_temp(n,:)*(1+2*A(n)^gamma(n)));
            R(n,:) = R(n,:) + temp(1:length(time));
        end
    end
    %% Convert spikes into activation
    [c,cf,A_tilde,A] = spike2activation(R(:,t),c,cf,A,parameter_Matrix,L_ce,S_i,Y_i,Fs);
    
    c_mat(:,t) = c;
    cf_mat(:,t) = cf;
    A_tilde_mat(:,t) = A_tilde;
    A_mat(:,t) = A;
    
    [F_ce(t),F_se(t)] = contraction_dynamics_v2(A,L_se,L_ce,V_ce,FL,FV,index_slow,Lmax,PTi,F0);
    
    k_0_de = h*MuscleVelocity(t);
    l_0_de = h*contraction_dynamics(A,L_se,L_ce,V_ce,FL,FV,modelParameter,index_slow,Lmax,PTi,F0);
    k_1_de = h*(MuscleVelocity(t)+l_0_de/2);
    l_1_de = h*contraction_dynamics(A,(Lmt - (L_ce+k_0_de/L0)*L0*cos(alpha))/L0T,L_ce+k_0_de/L0,V_ce+l_0_de/L0,FL,FV,modelParameter,index_slow,Lmax,PTi,F0);
    k_2_de = h*(MuscleVelocity(t)+l_1_de/2);
    l_2_de = h*contraction_dynamics(A,(Lmt - (L_ce+k_1_de/L0)*L0*cos(alpha))/L0T,L_ce+k_1_de/L0,V_ce+l_1_de/L0,FL,FV,modelParameter,index_slow,Lmax,PTi,F0);
    k_3_de = h*(MuscleVelocity(t)+l_2_de);
    l_3_de = h*contraction_dynamics(A,(Lmt - (L_ce+k_2_de/L0)*L0*cos(alpha))/L0T,L_ce+k_2_de/L0,V_ce+l_2_de/L0,FL,FV,modelParameter,index_slow,Lmax,PTi,F0);
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
output.v_mat = v_mat;
output.I_mat = I_mat;
output.Force = F_ce;
output.force = force;
output.ForceTendon = F_se;
output.Lce = MuscleLength./(L0/100);
output.Vce = MuscleVelocity./(L0/100);
output.Ace = MuscleAcceleration./(L0/100);
output.FR_Ia = FR_Ia;
output.FR_Ib = FR_Ib;
output.FR_RI = FR_RI;
output.U = U_mat;


%% Motoneuron
    function [u,v,spike_vec] = motoneuron(I,u,v,spike_vec,a,Fs)
        b_MN = 0.2;
        c_MN = -65;
        d_MN = 6;
        
        alpha_MN = 0.04;
        beta_MN  = 5;
        gamma_MN = 140;
        
        index_MN = find(v>=30);
        spike_vec(index_MN) = 1;
        v(index_MN) = c_MN;
        u(index_MN) = u(index_MN) + d_MN;
        
        v_dot = (alpha_MN*v.^2+beta_MN.*v+gamma_MN-u+I);
        v = v_dot*1000/Fs + v;
        
        u_dot = a.*(b_MN.*v-u);
        u = u_dot*1000/Fs + u;
        
    end
%% Convert spike trian into activation
    function [c,cf,A_tilde,A] = spike2activation(R,c,cf,A,parameter_Matrix,Lce,S_i,Y_i,Fs)
        S = parameter_Matrix(:,1);
        C = parameter_Matrix(:,2);
        k_1 = parameter_Matrix(:,3);
        k_2 = parameter_Matrix(:,4);
        k_3 = parameter_Matrix(:,5)*Lce + parameter_Matrix(:,6);
        k_4 = parameter_Matrix(:,7)*Lce + parameter_Matrix(:,8);
        tau_2 = parameter_Matrix(:,10);
        N = parameter_Matrix(:,11)*Lce + parameter_Matrix(:,12);
        K = parameter_Matrix(:,13)*Lce + parameter_Matrix(:,14);
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

%% Muscle Spindle
% Bag 1
    function [AP_bag1,f_dynamic,T,T_dot] = bag1_model(f_dynamic,gamma_dynamic,T,T_dot,L,V,A,Fs)
        p_bag_1 = 2;
        R_bag_1 = 0.46;
        a_bag_1 = 0.3;
        K_SR_bag_1 = 10.4649;
        K_PR_bag_1 = 0.15;
        M_bag_1 = 0.0002;
        LN_SR_bag_1 = 0.0423;
        L0_SR_bag_1 = 0.04;
        L0_PR_bag_1 = 0.76;
        tau_bag1 = 0.149;
        freq_bag1 = 60;
        beta0 = 0.0605;
        beta1 = 0.2592;
        Gamma1 = 0.0289;
        G = 20000;
        if V >= 0
            C = 1;
        else
            C = 0.42;
        end
        df_dynamic = (gamma_dynamic^p_bag_1/(gamma_dynamic^p_bag_1+freq_bag1^p_bag_1)-f_dynamic)/tau_bag1;
        f_dynamic = df_dynamic*1/Fs + f_dynamic;
        
        beta = beta0 + beta1 * f_dynamic;
        Gamma = Gamma1 * f_dynamic;
        
        T_ddot =  K_SR_bag_1/M_bag_1 * (C * beta * sign(V-T_dot/K_SR_bag_1)*((abs(V-T_dot/K_SR_bag_1))^a_bag_1)...
            *(L-L0_SR_bag_1-T/K_SR_bag_1-R_bag_1)+K_PR_bag_1*(L-L0_SR_bag_1-T/K_SR_bag_1-L0_PR_bag_1)+M_bag_1*A+Gamma-T);
        T_dot = T_ddot*1/Fs + T_dot;
        T = T_dot*1/Fs + T;
        
        AP_bag1 = G*(T/K_SR_bag_1-(LN_SR_bag_1-L0_SR_bag_1));
    end

%  Bag 2
    function [AP_primary_bag2,AP_secondary_bag2,f_static,T,T_dot] = bag2_model(f_static,gamma_static,T,T_dot,L,V,A,Fs)
        p = 2;
        R_bag_2 = 0.46;
        a_bag_2 = 0.3;
        K_SR_bag_2 = 10.4649;
        K_PR_bag_2 = 0.15;
        M_bag_2 = 0.0002;
        LN_SR_bag_2 = 0.0423;
        LN_PR_bag_2 = 0.89;
        L0_SR_bag_2 = 0.04;
        L0_PR_bag_2 = 0.76;
        L_secondary = 0.04;
        X = 0.7;
        tau_bag2 = 0.205;
        freq_bag2 = 60;
        beta0 = 0.0822;
        beta2 = -0.046;
        Gamma2 = 0.0636;
        G = 10000;
        if V >= 0
            C = 1;
        else
            C = 0.42;
        end
        df_static = (gamma_static^p/(gamma_static^p+freq_bag2^p)-f_static)/tau_bag2;
        f_static = df_static*1/Fs + f_static;
        beta = beta0 + beta2 * f_static;
        Gamma = Gamma2 * f_static;
        
        T_ddot = K_SR_bag_2/M_bag_2 * (C * beta * sign(V-T_dot/K_SR_bag_2)*((abs(V-T_dot/K_SR_bag_2))^a_bag_2)...
            *(L-L0_SR_bag_2-T/K_SR_bag_2-R_bag_2)+K_PR_bag_2*(L-L0_SR_bag_2-T/K_SR_bag_2-L0_PR_bag_2)+M_bag_2*A+Gamma-T);
        T_dot = T_ddot*1/Fs + T_dot;
        T = T_dot*1/Fs + T;
        
        AP_primary_bag2 = G*(T/K_SR_bag_2-(LN_SR_bag_2-L0_SR_bag_2));
        AP_secondary_bag2 = G*(X*L_secondary/L0_SR_bag_2*(T/K_SR_bag_2-(LN_SR_bag_2-L0_SR_bag_2))+(1-X)*L_secondary/L0_PR_bag_2*(L-T/K_SR_bag_2-L0_SR_bag_2-LN_PR_bag_2));
        
    end

% Chain
    function [AP_primary_chain,AP_secondary_chain,T,T_dot] = chain_model(gamma_static,T,T_dot,L,V,A,Fs)
        p = 2;
        R_sp = 0.46;
        a = 0.3;
        K_SR = 10.4649;
        K_PR= 0.15;
        M = 0.0002;
        LN_SR = 0.0423;
        LN_PR = 0.89;
        L0_SR = 0.04;
        L0_PR = 0.76;
        L_secondary = 0.04;
        X = 0.7;
        freq_chain = 90;
        beta0 = 0.0822;
        beta2_chain = -0.069;
        Gamma2_chain = 0.0954;
        G_chain = 10000;
        if V >= 0
            C = 1;
        else
            C = 0.42;
        end
        
        f_static_chain = gamma_static^p/(gamma_static^p+freq_chain^p);
        beta_chain = beta0 + beta2_chain * f_static_chain;
        Gamma_chain = Gamma2_chain * f_static_chain;
        
        T_ddot_chain = K_SR/M * (C * beta_chain * sign(V-T_dot/K_SR)*((abs(V-T_dot/K_SR))^a)...
            *(L-L0_SR-T/K_SR-R_sp)+K_PR*(L-L0_SR-T/K_SR-L0_PR)+M*A+Gamma_chain-T);
        T_dot = T_ddot_chain*1/Fs + T_dot;
        T = T_dot*1/Fs + T;
        
        AP_primary_chain = G_chain*(T/K_SR-(LN_SR-L0_SR));
        AP_secondary_chain = G_chain*(X*L_secondary/L0_SR*(T/K_SR-(LN_SR-L0_SR))+(1-X)*L_secondary/L0_PR*(L-T/K_SR-L0_SR-LN_PR));
        
    end

% Generate primary and secondary output
    function [Output_Primary,Output_Secondary] = SpindleOutput(AP_bag1,AP_primary_bag2,AP_secondary_bag2,AP_primary_chain,AP_secondary_chain)
        S = 0.156;
        
        if AP_bag1 < 0
            AP_bag1 = 0;
        end
        
        if AP_primary_bag2 < 0
            AP_primary_bag2 = 0;
        end
        
        if AP_primary_chain < 0
            AP_primary_chain = 0;
        end
        
        
        if AP_secondary_bag2 < 0
            AP_secondary_bag2 = 0;
        end
        
        if AP_secondary_chain < 0
            AP_secondary_chain = 0;
        end
        
        
        if AP_bag1 > (AP_primary_bag2+AP_primary_chain)
            Larger = AP_bag1;
            Smaller = AP_primary_bag2+AP_primary_chain;
        elseif AP_bag1 < (AP_primary_bag2+AP_primary_chain)
            Larger = AP_primary_bag2+AP_primary_chain;
            Smaller = AP_bag1;
        elseif AP_bag1 == (AP_primary_bag2+AP_primary_chain)
            Larger = 0;
            Smaller = 0;
        end
        
        Output_Primary = Larger + S * Smaller;
        Output_Secondary = AP_secondary_bag2 + AP_secondary_chain;
        
        if Output_Primary < 0
            Output_Primary = 0;
        elseif Output_Primary > 100000
            Output_Primary = 100000;
        end
        if Output_Secondary < 0
            Output_Secondary = 0;
        elseif Output_Secondary > 100000
            Output_Secondary = 100000;
        end
        
    end

    function ddx = contraction_dynamics(A,L_s,L_m,L_m_dot,FL_vec,FV_vec,modelParameter,index_slow,Lmax,PT,F0)
        %% Force-length and force-velocity
        FL_vec(1:index_slow) = FL_slow_function(L_m);
        FL_vec(index_slow+1:end) = FL_fast_function(L_m);
        
        if L_m_dot > 0
            FV_vec(1:index_slow) = FVecc_slow_function(L_m,L_m_dot);
            FV_vec(index_slow+1:end) = FVecc_fast_function(L_m,L_m_dot);
        else
            FV_vec(1:index_slow) = FVcon_slow_function(L_m,L_m_dot);
            FV_vec(index_slow+1:end) = FVcon_fast_function(L_m,L_m_dot);
        end
        
        %% Passive element 1
        F_pe1 = Fpe1_function(L_m/Lmax,L_m_dot);
        
        %% Passive element 2
        F_pe2 = Fpe2_function(L_m);
        if F_pe2 > 0
            F_pe2 = 0;
        end
        
        f_i = A.*PT.*(FL_vec.*FV_vec+F_pe2);
        
        F_m_temp = sum(f_i);
        F_m = F_m_temp + F_pe1*F0;
        
        F_t = Fse_function(L_s) * F0;
        
        M = modelParameter.mass;
        rho = modelParameter.pennationAngle;
        
        ddx = (F_t*cos(rho) - F_m*(cos(rho)).^2)/(M) ...
            + (L_m_dot).^2*tan(rho).^2/(L_m);
    end

    function [F_m,F_t] = contraction_dynamics_v2(A,L_s,L_m,L_m_dot,FL_vec,FV_vec,index_slow,Lmax,PT,F0)
        %% Force-length and force-velocity
        FL_vec(1:index_slow) = FL_slow_function(L_m);
        FL_vec(index_slow+1:end) = FL_fast_function(L_m);
        
        if L_m_dot > 0
            FV_vec(1:index_slow) = FVecc_slow_function(L_m,L_m_dot);
            FV_vec(index_slow+1:end) = FVecc_fast_function(L_m,L_m_dot);
        else
            FV_vec(1:index_slow) = FVcon_slow_function(L_m,L_m_dot);
            FV_vec(index_slow+1:end) = FVcon_fast_function(L_m,L_m_dot);
        end
        
        %% Passive element 1
        F_pe1 = Fpe1_function(L_m/Lmax,L_m_dot);
        
        %% Passive element 2
        F_pe2 = Fpe2_function(L_m);
        if F_pe2 > 0
            F_pe2 = 0;
        end
        
        f_i = A.*PT.*(FL_vec.*FV_vec+F_pe2);
        
        F_m_temp = sum(f_i);
        F_m = F_m_temp + F_pe1*F0;
        
        F_t = Fse_function(L_s) * F0;
        
    end

    function [x] = noise(x,D,Fs)
        vec_length = size(x,1);
        tau = 0.01;
        chi = normrnd(0,1,[vec_length,1]);
        x_dot = -x./tau + sqrt(D)*chi;
        x = x_dot.*1/Fs + x;
    end

end