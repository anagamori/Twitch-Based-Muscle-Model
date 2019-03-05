%==========================================================================
% new_model_muscle_parameters.m
% Author: Akira Nagamori
% Last update: 3/4/19
% Descriptions:
%   Define parameters related to muscle 
%       - muscle architecture (e.g., muscle length, tendon length)
%       - tetanic tension of whole muscle and individual motor units
%       - distribution of contraction time across motor units
%==========================================================================

%% Muscle architectural parameters
modelParameter.pennationAngle = 9.6*pi/180; %[radians]
modelParameter.optimalLength = 6.8; % [cm]
modelParameter.tendonSlackLength = 24.1; % [cm]
modelParameter.mass = 0.01; % [g]
modelParameter.muscleInitialLength = 6.8; % [cm]
modelParameter.tendonInitialLength = 24.1; % [cm]

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

%--------------------------------------------------------------------------
% Motor unit architecture
N_MU = 300; % number of motor units
i_MU = 1:N_MU; % index for motor units

%% Peak tetanic force
RP_MU = 25; %range of peak tension across motor untis in unit of fold
b_MU = log(RP_MU)/N_MU; %coefficient to establish a range of twich force values
P_MU = exp(b_MU*i_MU); %force generated by a motor unit as a function of its recruitment threshold
PTi = P_MU./sum(P_MU)*F0; % peak tetanic force for individual units

%% Fractional PSCA
F_pcsa_slow = 0.3; % fractional PSCA of slow-twitch motor units (0-1)
[~, index_slow] = min(abs(cumsum(PTi) - F0*F_pcsa_slow)); rng(1)

%% Contraction time
% Generate a distribution of contraction time across motor units based on
% Rayleigh distribution
rng(1)
min_CT = 20; %minimum contraction time [ms]
CT = round(raylrnd(23,1,N_MU)+min_CT); %contraction time of individual motor units [ms]
histogram(CT)
CT_sorted = sort(CT,'descend');
CT_fastest_slow = CT_sorted(index_slow);
mean(CT_sorted(1:index_slow)); %average contraction time of slow-twitch MUs
mean(CT_sorted(index_slow+1:end)); %average contraction time of fast-twitch MUs

load('CT_vec')
[CT_new,index_CT] = sort(CT_vec,'descend');
%% Motor unit parameter
load('pool_parameter_matrix')
parameter_Matrix; %matrix of [N_MU,15]
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
% plot(CT_sorted,PTi(index_MU_PTi),'o')

%% Recruitment threshold
% Find recruitment threshold for individual units using exponential fit
% Recruitment threshold is correlated to peak tetanic tension 
%   Use index_MU_PTi to appropriately index each MU
Ur = 0.8; % recruitment threshold for the lastly recruited motor unit
Ur_1 = 0.01; % reruitment threshold for the first unit
f_RT = fit([1 N_MU]',[Ur_1 Ur]','exp1');
coeffs_f_RT = coeffvalues(f_RT);
U_th = coeffs_f_RT(1)*exp(coeffs_f_RT(2)*i_MU); % the resulting recruitment threshold for individual units

%% Minimum and maximum firing rate
load('FR_half')
FR_half_new = FR_half(index_CT);
MFR = FR_half/2;
PFR = FR_half*2;

%%

