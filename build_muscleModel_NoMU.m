%==========================================================================
% build_muscleModel_NoMU.m
% Author: Akira Nagamori
% Last update: 3/25/19
% Descriptions:
%   Create a model with 100 motor units in a pool
%   Maintain the maximum force of muscle, which force individual units to
%   produce larger peak tetanic force
%   Sample 100 units from 300 units in the base model 
%   Save model parameters as data structure (modelParameter) into an appropriate folder
%==========================================================================
close all
clear all
clc

code_folder = '/Users/akira/Documents/Github/Twitch-Based-Muscle-Model';
model_parameter_folder =  '/Users/akira/Documents/Github/Twitch-Based-Muscle-Model/Model Parameters/Model_N_100_sameF0';
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
modelParameter.F0 = PCSA * sigma; % maximal force

L_tendon = modelParameter.tendonSlackLength;
modelParameter.L0T = L_tendon*1.05;
alpha = modelParameter.pennationAngle;
Lm_initial = modelParameter.muscleInitialLength; % muscle initial length
Lt_initial = modelParameter.tendonInitialLength; % tendon initial length
modelParameter.Lmt = Lm_initial*cos(alpha)+Lt_initial; % intial musculotendon length
[modelParameter.L_ce,modelParameter.L_se,modelParameter.Lmax] =  InitialLength_function(modelParameter);

%% Motor unit parameters
modelParameter.N_MU = 100; % number of motor units
modelParameter.i_MU = 1:modelParameter.N_MU; % index for motor units

%% Peak tetanic force
RP_MU = 25; %range of peak tension across motor untis in unit of fold
b_MU = log(RP_MU)/modelParameter.N_MU; %coefficient to establish a range of twich force values
P_MU = exp(b_MU*modelParameter.i_MU); %force generated by a motor unit as a function of its recruitment threshold
PTi = P_MU./sum(P_MU)*modelParameter.F0; % peak tetanic force for individual units

%% Fractional PSCA
F_pcsa_slow = 0.3; % fractional PSCA of slow-twitch motor units (0-1)
[~, modelParameter.index_slow] = min(abs(cumsum(PTi) - modelParameter.F0*F_pcsa_slow)); 

%% Model parameters for activation-frequency relationship
cd(model_parameter_folder )
load('pool_parameter_matrix')
modelParameter.parameterMatrix = parameter_Matrix(1:300/modelParameter.N_MU:300,:);
cd(code_folder)

%% Assign peak tetanic force into each unit
% shuffle indexes within each fiber type with respect to contraction time
% this will allow us to randomly assign peak tetanic tension to each motor
% unit with different contraction time
rng(1)
R_slow = randperm(modelParameter.index_slow);
index_fast = modelParameter.index_slow+1:modelParameter.N_MU;
R_fast_temp = randperm(length(index_fast));
R_fast = index_fast(R_fast_temp);
index_MU_PTi = [R_slow R_fast]; % vector of indexes to match peak tetanic tension to appropriate contraction time
modelParameter.PTi_new = PTi (index_MU_PTi);

%% Recruitment threshold
% Find recruitment threshold for individual units using exponential fit
% Recruitment threshold is correlated to peak tetanic tension
%   Use index_MU_PTi to appropriately index each MU
Ur = 0.5; % recruitment threshold for the lastly recruited motor unit
Ur_1 = 0.01; % reruitment threshold for the first unit
f_RT = fit([1 modelParameter.N_MU]',[Ur_1 Ur]','exp1');
coeffs_f_RT = coeffvalues(f_RT);
U_th = coeffs_f_RT(1)*exp(coeffs_f_RT(2)*modelParameter.i_MU); % the resulting recruitment threshold for individual units
modelParameter.U_th_new = U_th(index_MU_PTi);

%% Minimum and maximum firing rate
cd(model_parameter_folder )
load('FR_half')
cd(code_folder)
modelParameter.FR_half = FR_half(1:300/modelParameter.N_MU:300);
modelParameter.MDR = modelParameter.FR_half/2;
modelParameter.PDR = modelParameter.FR_half*2;

%% Gain for frequency-activation relationship
modelParameter.g_e = max((modelParameter.PDR-modelParameter.MDR)./(1-modelParameter.U_th_new));

%% Sample 100 units from the base model

%% Save model parameters
cd(model_parameter_folder)
save('modelParameter','modelParameter')
cd(code_folder)

%%
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
