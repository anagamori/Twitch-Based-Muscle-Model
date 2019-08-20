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
Ia_Input = zeros(1,length(time));

%%
Ib_gain = SLRParameter.Ib_gain;
Ib_delay = SLRParameter.Ib_delay;

x_GTO = zeros(1,length(time));
FR_Ib_temp = zeros(1,length(time));
FR_Ib = zeros(1,length(time));
Ib_Input = zeros(1,length(time));

%%
RI_gain = SLRParameter.Ib_gain;
RI_delay = SLRParameter.RI_delay;

FR_RI_temp = zeros(1,length(time));
FR_RI = zeros(1,length(time));
RI_Input = zeros(1,length(time));

%% 
C_delay = SLRParameter.C_delay;
K_C = SLRParameter.K_C;
%% 
noise_Ia = 0;
noise_Ib = 0;
noise_RI = 0;
noise_C = 0;
noise_ID = zeros(N_MU,1);
noise_CM = zeros(1,1);
%%
C_temp = 0;
U_vec = zeros(1,length(time));
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
    Ia_Input(t) = FR_Ia(t)/Ia_gain;
    [noise_Ia] = noise(noise_Ia,0,Fs);
    Ia_Input(t) = Ia_Input(t) + Ia_Input(t)*noise_Ia;
    %%
    if t > 5
        [FR_Ib,FR_Ib_temp,x_GTO] = GTOOutput(FR_Ib,FR_Ib_temp,x_GTO,F_se(t-1),t);
        [FR_RI,FR_RI_temp] = RenshawOutput(FR_RI,FR_RI_temp,U_eff,t);
    end
    Ib_Input(t) = FR_Ib(t)/Ib_gain;
    [noise_Ib] = noise(noise_Ib,0,Fs);
    Ib_Input(t) = Ib_Input(t) + Ib_Input(t)*noise_Ib;
    
    RI_Input(t) = FR_RI(t)/RI_gain;
    [noise_RI] = noise(noise_RI,0,Fs);
    RI_Input(t) = RI_Input(t) + RI_Input(t)*noise_RI;
    
    %%
    [noise_C] = noise(noise_C,10000,Fs);
    
    [noise_ID] = noise(noise_ID,10000,Fs);
    [noise_CM] = noise(noise_CM,0,Fs);
    %%
    if t > 1
        if controlOpt == 1
           U = C_input(t);
        elseif controlOpt == 2
            %% Effective activation (Song et al., 2008)
            if t > Ia_delay && t <= Ib_delay
                U = C_input(t) + noise_C*C_input(t) ...
                    + Ia_Input(t-Ia_delay) ...
                    - RI_Input(t-RI_delay);
            elseif t > Ib_delay && t <= C_delay 
                U = C_input(t) + noise_C*C_input(t) ...
                    + Ia_Input(t-Ia_delay) ...
                    - Ib_Input(t-Ib_delay)...
                    - RI_Input(t-RI_delay);
            elseif t > C_delay 
                C_temp = K_C*(C_input(t) - F_se(t-C_delay)/F0) + C_temp; 
                U = C_temp + noise_C*C_temp...
                    + Ia_Input(t-Ia_delay) ...
                    - Ib_Input(t-Ib_delay)...
                    - RI_Input(t-RI_delay);
            else
                U = C_input(t) + noise_C*C_input(t);
            end
            if U < 0 
                U = 0;
            end
        end
        U = U+noise_CM*U;      
        U_vec(t) = U; 
        
        %% Calculate firing rate
        % Linear increase in discharge rate up to Ur
        if recruitmentType == 1
            I = g_e.*(U+noise_ID*U-U_th) + I_th;
        elseif recruitmentType == 3
            I = zeros(N_MU,1);
            U_temp = U + noise_ID*U;
            I_temp_1 = I_th + lamda.*k_e.*(U_temp-U_th);
            index_1 = find(U_temp <= U_th_t);
            I(index_1) = I_temp_1(index_1);
            I_temp_2 = I_max-k_e.*(1-U_temp);
            index_2 = find(U_temp > U_th_t);
            I(index_2) = I_temp_2(index_2);
            I_temp_3 = g_e.*(U_temp-U_th)+I_th;
            I(index_t) = I_temp_3(index_t);
            %I = (I+I.*noise_ID'+I*noise_CM)';
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
    
    f_i = A.*PTi'.*(FL.*FV+F_pe2);
    %f_i = A.*PTi_new'.*(FL+F_pe2);
    force(:,t) = f_i;
    
    F_ce(t) = sum(f_i);
    F_total(t) = F_ce(t) + F_pe1*F0;
    
    F_se(t) = Fse_function(L_se) * F0;
    
    MuscleAcceleration(t+1) =(F_se(t)*cos(alpha) - F_total(t)*(cos(alpha)).^2)/(mass) ...
        + (MuscleVelocity(t)).^2*tan(alpha).^2/(MuscleLength(t));
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
    A_ce = MuscleAcceleration(t+1)/(L0/100);
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
output.U = U_vec;


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

    %%  GTO
    function [FR_Ib,FR_Ib_temp,x_GTO] = GTOOutput(FR_Ib,FR_Ib_temp,x_GTO,Force,index)
        G1 = 60;
        G2 = 4;
        num1 = 1.7;
        num2 = -3.399742022978487;
        num3 = 1.699742026978047;
        den1 = 1.0;
        den2 = -1.999780020198665;
        den3 = 0.999780024198225;
                
        x_GTO(index) = G1*log(Force/G2+1);
              
        FR_Ib_temp(index) = (num3*x_GTO(index-2) + num2*x_GTO(index-1) + num1*x_GTO(index) - ...
                  den3*FR_Ib_temp(index-2) - den2*FR_Ib_temp(index-1))/den1;
        FR_Ib(index) = FR_Ib_temp(index);
        if FR_Ib(index)<0
            FR_Ib(index) = 0;
        end
    end
    
    %% Renshaw cell
    function [FR_RI,FR_RI_temp] = RenshawOutput(FR_RI,FR_RI_temp,ND,index)
        num1 = 0.238563173450928;
        num2 = -0.035326319453965;
        num3 = -0.200104635331441;
        den1 = 1.0;
        den2 = -1.705481699867712;
        den3 = 0.708613918533233;
        
        FR_RI_temp(index) = (num3*ND(index-2)+num2*ND(index-1)+num1*ND(index)-den3*FR_RI_temp(index-2)-den2*FR_RI_temp(index-1))/den1;
        FR_RI(index) = FR_RI_temp(index); 
        if FR_RI(index) < 0
            FR_RI(index) = 0;
        end
        
    end

    function [x] = noise(x,D,Fs)
        vec_length = size(x,1);
        tau = 0.01;
        chi = normrnd(0,1,[vec_length,1]);
        x_dot = -x./tau + sqrt(D)*chi;
        x = x_dot.*1/Fs + x;
    end

end