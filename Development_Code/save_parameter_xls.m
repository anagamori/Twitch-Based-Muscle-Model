close all
clc
clear all

code_folder = '/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Development_Code/';
model_parameter_folder =  '/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Development_Code/Data';
cd(model_parameter_folder)
load('modelParameter_v2')
cd(code_folder)

parameterMatrix = modelParameter.parameterMatrix;

MUID = 1:modelParameter.N_MU;
MUID = MUID';
S = round(parameterMatrix(:,1),2);
C = round(parameterMatrix(:,2),2);
k_1 = round(parameterMatrix(:,3),2);
k_2 = round(parameterMatrix(:,4),2);
k_3 = round(parameterMatrix(:,5),2);
k_4 = round(parameterMatrix(:,6),2);
tau_1 = round(parameterMatrix(:,7),4);
tau_2 = round(parameterMatrix(:,8),4);
N = round(parameterMatrix(:,9),2);
K = round(parameterMatrix(:,10),4);
tau_3 = round(parameterMatrix(:,11),4);
gamma = round(parameterMatrix(:,12),2);
phi_1 = round(parameterMatrix(:,13),2);
phi_2 = round(parameterMatrix(:,14),2);

T = table(MUID,S,C,k_1,k_2,k_3,k_4,tau_1,tau_2,N,K,tau_3,gamma,phi_1,phi_2);

filename = 'MU_parameters.xlsx';
cd(model_parameter_folder)
writetable(T,filename,'Sheet',1)
cd(code_folder)