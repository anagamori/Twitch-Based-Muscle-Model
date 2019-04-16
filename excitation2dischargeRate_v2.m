%==========================================================================
% excitation2dischargeRate.m
% Author: Akira Nagamori
% Last update: 3/5/19
% Descriptions:
%   Plot the relationship between excitation and discharge rate of
%   individual motor units
%==========================================================================

cd('/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Model Parameters/Model_CTvsPTi');
%% Peak tension of muscle
density = 1.06; %
L0 = 6.8; % optimal muscle length [cm]
mass = 0.01; % muscle mass [kg]
PCSA = (mass*1000)/(density*L0); % PCSA of muscle
sigma = 31.8; % specific tension
F0 = PCSA * sigma; % maximal force

%% Number of motor unit
N_MU = 300; % number of motor units
i_MU = 1:N_MU; % index for motor units

%% Contraction time
% Generate a distribution of contraction time across motor units based on
% Rayleigh distribution
% rng(1)
% min_CT = 20; %minimum contraction time [ms]
% CT = round(raylrnd(23,1,N_MU)+min_CT); %contraction time of individual motor units [ms]
% CT_sorted = sort(CT,'descend');
load('CT_vec')
modelParameter.CT = CT_vec;

%% Peak tetanic force
RP_MU = 25; %range of peak tension across motor untis in unit of fold
b_MU = log(RP_MU)/N_MU; %coefficient to establish a range of twich force values
P_MU = exp(b_MU*i_MU); %force generated by a motor unit as a function of its recruitment threshold
PTi = P_MU./sum(P_MU)*F0; % peak tetanic force for individual units

%% Fractional PSCA
F_pcsa_slow = 0.3; % fractional PSCA of slow-twitch motor units (0-1)
[~, index_slow] = min(abs(cumsum(PTi) - F0*F_pcsa_slow)); rng(1)
%index_slow = 196;

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
PTi_new = PTi; %(index_MU_PTi);

%% Recruitment threshold
% Find recruitment threshold for individual units using exponential fit
% Recruitment threshold is correlated to peak tetanic tension
%   Use index_MU_PTi to appropriately index each MU
Ur = 0.8; % recruitment threshold for the lastly recruited motor unit
Ur_1 = 0.01; % reruitment threshold for the first unit
f_RT = fit([1 N_MU]',[Ur_1 Ur]','exp1');
coeffs_f_RT = coeffvalues(f_RT);
U_th = coeffs_f_RT(1)*exp(coeffs_f_RT(2)*i_MU); % the resulting recruitment threshold for individual units
U_th_new = U_th; %(index_MU_PTi);
[~,loc_max_U_th] = max(U_th_new);

%% FR_half for individual motor units
load('FR_half')
MDR = FR_half/2;
PDR = FR_half*2;

[~,index_DR_dif] = max(PDR-MDR);
%g_e = (PDR(index_DR_dif)-MDR(index_DR_dif))/(1-U_th_new(index_DR_dif));
g_e = (2-0.5)./(1-U_th_new(end));
%g_e = 115.1750;
%% 
cd('/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model')
%% Discharge rate of motor unit
U_vec = 0:0.01:1;
DR_mat = zeros(N_MU,length(U_vec));
for i = 1:length(U_vec)
    DR_temp = g_e.*(U_vec(i)-U_th_new)+0.5;
    DR_temp(DR_temp<0.5) = 0;
    DR_temp(DR_temp>2) = 2;
    
    DR_MU = DR_temp.*FR_half;
    DR_mat(:,i) = DR_MU;
end

figure(1)
plot(U_vec,DR_mat)
xlabel('Synaptic Drive (%Maximum)')
ylabel('Discharge Rate (Hz)')

