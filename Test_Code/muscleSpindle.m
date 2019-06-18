function [outputIa,outputII,output_bag1,output_bag2,output_chain] = muscleSpindle(Lce,Vce,Ace,Fs,gamma_dynamic,gamma_static)
%%%%%%%%%%% Muscle Spindle Model
% Input: Lce = muscle legnth normalized to optimal muscle legnth (L0)
%        Vce = rate of change in muscle length, or velcoity, in unit of L0/s
%        Ace = acceleration of muscle length in unit of L0/s^2
%        Fs = sampling frequency 
%        gamma_dynamic = dynamic gamma drive in unit of pps
%        gamma_static = static gamma drive in unit of pps
% Output: Ia afferent firing rate = outputIa
%         II afferent firing rate = outputII

% Time 
t = 0:1/Fs:(length(Lce)-1)/Fs;

%%%%%%%%%%% Spindle Model Parameters %%%%%%%%%%%%%%
% all the parameters can be found in Mileusnic et al. (2006)

p = 2;
R = 0.46; 
a = 0.3;
K_SR = 10.4649;
K_PR = 0.15;
M = 0.0002; 

LN_SR = 0.0423;
LN_PR = 0.89;

L0_SR = 0.04;
L0_PR = 0.76; 

L_secondary = 0.04;
X = 0.7;

tau_bag1 = 0.149;
freq_bag1 = 60;

beta0_bag1 = 0.0605;
beta1 = 0.2592;
Gamma1 = 0.0289;

G_bag1 = 20000; 
tau_bag2 = 0.205;
freq_bag2 = 60;

beta0 = 0.0822;
beta2 = -0.046;
Gamma2 = 0.0636;

G_bag2 = 10000;

freq_chain = 90;

beta0 = 0.0822;
beta2_chain = - 0.069;
Gamma2_chain = 0.0954;

G_chain = 10000;

f_dynamic = 0;
f_static = 0;
T_ddot_bag1 = 0;
T_dot_bag1 = 0;
T_bag1 = 0;
T_ddot_bag2 = 0;
T_dot_bag2 = 0;
T_bag2 = 0;
T_ddot_chain = 0;
T_dot_chain = 0;
T_chain = 0;

%%%%%%%%% Parameter Initialization 
outputIa = zeros(1,length(Lce));
outputII = zeros(1,length(Lce));

%%%%%%%%


for i = 1:length(Lce)
    L = Lce(i);
    L_dot = Vce(i);
    L_ddot = Ace(i);
    
    
    if Vce(i) >= 0
        C = 1;
    else
        C = 0.42;
    end
    
    %%
    df_dynamic = (gamma_dynamic^p/(gamma_dynamic^p+freq_bag1^p)-f_dynamic)/tau_bag1;
    f_dynamic = 1/Fs*df_dynamic + f_dynamic;
    
    beta_bag1 = beta0_bag1 + beta1 * f_dynamic;
    Gamma_bag1 = Gamma1 * f_dynamic;
    
    T_ddot_bag1 = K_SR/M * (C * beta_bag1 * sign(L_dot-T_dot_bag1/K_SR)*((abs(L_dot-T_dot_bag1/K_SR))^a)...
        *(L-L0_SR-T_bag1/K_SR-R)+K_PR*(L-L0_SR-T_bag1/K_SR-L0_PR)+M*L_ddot+Gamma_bag1-T_bag1);
    T_dot_bag1 = T_ddot_bag1*1/Fs + T_dot_bag1;
    T_bag1 = T_dot_bag1*1/Fs + T_bag1;
    
    AP_bag1 = G_bag1*(T_bag1/K_SR-(LN_SR-L0_SR));
    
    %%
    df_static = (gamma_static^p/(gamma_static^p+freq_bag2^p)-f_static)/tau_bag2;
    f_static = 1/Fs*df_static + f_static;
    
    beta_bag2 = beta0 + beta2 * f_static;
    Gamma_bag2 = Gamma2 * f_static;
    
    T_ddot_bag2 = K_SR/M * (C * beta_bag2 * sign(L_dot-T_dot_bag2/K_SR)*((abs(L_dot-T_dot_bag2/K_SR))^a)...
        *(L-L0_SR-T_bag2/K_SR-R)+K_PR*(L-L0_SR-T_bag2/K_SR-L0_PR)+M*L_ddot+Gamma_bag2-T_bag2);
    T_dot_bag2 = T_ddot_bag2*1/Fs + T_dot_bag2;
    T_bag2 = T_dot_bag2*1/Fs + T_bag2;
    
    AP_primary_bag2 = G_bag2*(T_bag2/K_SR-(LN_SR-L0_SR));
    AP_secondary_bag2 = G_bag2*(X*L_secondary/L0_SR*(T_bag2/K_SR-(LN_SR-L0_SR))+(1-X)*L_secondary/L0_PR*(L-T_bag2/K_SR-L0_SR-LN_PR));
    
    %%
    f_static_chain = gamma_static^p/(gamma_static^p+freq_chain^p);
    
    beta_chain = beta0 + beta2_chain * f_static_chain;
    Gamma_chain = Gamma2_chain * f_static_chain;
    
    T_ddot_chain = K_SR/M * (C * beta_chain * sign(L_dot-T_dot_chain/K_SR)*((abs(L_dot-T_dot_chain/K_SR))^a)...
        *(L-L0_SR-T_chain/K_SR-R)+K_PR*(L-L0_SR-T_chain/K_SR-L0_PR)+M*L_ddot+Gamma_chain-T_chain);
    T_dot_chain = T_ddot_chain*1/Fs + T_dot_chain;
    T_chain = T_dot_chain*1/Fs + T_chain;
    
    AP_primary_chain = G_chain*(T_chain/K_SR-(LN_SR-L0_SR));
    AP_secondary_chain = G_chain*(X*L_secondary/L0_SR*(T_chain/K_SR-(LN_SR-L0_SR))+(1-X)*L_secondary/L0_PR*(L-T_chain/K_SR-L0_SR-LN_PR));
    
    
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
    OutputPrimary = Larger + S * Smaller;
    OutputSecondary = AP_secondary_bag2 + AP_secondary_chain;
    
    if OutputPrimary < 0
        OutputPrimary = 0;
    elseif OutputPrimary > 100000
        OutputPrimary = 100000;
    end
    if OutputSecondary < 0
        OutputSecondary = 0;
    elseif OutputSecondary > 100000
        OutputSecondary = 100000;
    end
    outputIa(i) = OutputPrimary; % Ia afferent firing rate
    outputII(i) = OutputSecondary; % II afferent firing rate
    output_bag1(i) = AP_bag1;
    output_bag2(i) = AP_primary_bag2; % II afferent firing rate
    output_chain(i) = AP_primary_chain; % II afferent firing rate
end


figure()
subplot(3,1,1)
plot(t,outputIa)
ylabel('Ia Afferent Firing Rate')
subplot(3,1,2)
plot(t,outputII)
ylabel('II Afferent Firing Rate')
subplot(3,1,3)
plot(t,Lce)
ylabel('Muscle Length')

end



