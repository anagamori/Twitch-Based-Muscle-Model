%==========================================================================
% spikeDriveMuscleModel_parameterFit_FT.m
% Author: Akira Nagamori
% Last update: 3/27/19
% Descriptions: 
%   Adjust model parameters to a unit with a target contraction time
%==========================================================================

close all
clear all
clc

%% Folder name
code_folder = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model';
data_folder = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Model Parameters/Model_1';

%% Simulation parameters
Fs = 1000; %sampling frequency
T = 1/Fs;
time = 0:1/Fs:5; %simulation time

Lce = 1; % Muscle length

%% Load a vector that contains contraction times of each unit
cd (data_folder)
load('CT')
cd (code_folder)

%% Load seed parameters
cd(data_folder)
load('Seed_FT_2')
cd(code_folder)
param_seed = param;
clear param
%%
first_MU = 260; %150;
last_MU = 260; %179;
Data_cell = cell(1,last_MU);

%% Weighting for optimization
rng shuffle
b = 30;
a = 20;
weight_temp = ((b-a)*rand(1,last_MU)+a);
%parpool(10)

%%
param_target = 1:1:11;
M = randomizer(param_target,6,last_MU);
%% Test each unit
%parpool(10)
for j = first_MU:last_MU
    j
    param = param_seed;
    %% Contraction time of the unit to be optimized to
    target_CT = CT_sorted(j);
    index_temp = M(j,:);
    %% Run annealing curve algorithm for 6 iterations
    tic
    for k = 1:6
        
        Param_matrix = annealing_curve(param,k);
       
        %% Loop through all parameters
        for n = 1:11
            index = index_temp(n+(k-1)*11);
            %% Loop through all perturbations
            error_long = zeros(1,3);
            param_temp = zeros(3,11);
            for i = 1:3
                [S,C,k_1,k_2,k_3,k_4,tau_1,tau_2,N,K,alpha] = parameter_Assigning(param,Param_matrix,index,i);
                param_temp(i,:) =  [S,C,k_1,k_2,k_3,k_4,tau_1,tau_2,N,K,alpha];
            end
            for l = 1:3          
                %% Run a twitch simulation and sweep simulation
                [Data_temp_2] = spikeDrivenMuscleModel_testFunction(param_temp(l,:),Lce,0,'fast',0);
                weight = weight_temp(j);
                error_long(l) = Data_temp_2{2,8} + abs(target_CT-Data_temp_2{2,1}) + weight*Data_temp_2{2,5};
                
            end
            [min_error,loc_min_error] = min(error_long);
            [S,C,k_1,k_2,k_3,k_4,tau_1,tau_2,N,K,alpha] = parameter_Assigning(param,Param_matrix,index,loc_min_error);
            param =  [S,C,k_1,k_2,k_3,k_4,tau_1,tau_2,N,K,alpha];
        end
    end
    toc
    [Data_temp] = spikeDrivenMuscleModel_testFunction(param,Lce,0,'fast',1);
    Data_cell{j} = Data_temp;
    
end

for MU_No = first_MU:last_MU
    Data = Data_cell{MU_No};
    cd('/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Model Parameters/Model_1/FT')
    save(['Data_' num2str(MU_No)],'Data')
    cd(code_folder)
    Data{2,5}
end

%%

function Y = annealing_curve(x,k)

perturbation_amp = 0.3/2;

Y(1,1) = x(1) +  (x(1)*perturbation_amp)./2.^(k-2);
Y(2,1) = x(1) -  (x(1)*perturbation_amp)./2.^(k-2);
Y(3,1) = x(1);

Y(1,2) = x(2) +  (x(2)*perturbation_amp)./2.^(k-2);
Y(2,2) = x(2) -  (x(2)*perturbation_amp)./2.^(k-2);
Y(3,2) = x(2);

Y(1,3) = x(3) +  (x(3)*perturbation_amp)./2.^(k-2);
Y(2,3) = x(3) -  (x(3)*perturbation_amp)./2.^(k-2);
Y(3,3) = x(3);

Y(1,4) = x(4) +  (x(4)*perturbation_amp)./2.^(k-2);
Y(2,4) = x(4) -  (x(4)*perturbation_amp)./2.^(k-2);
Y(3,4) = x(4);

Y(1,5) = x(5) +  (x(5)*perturbation_amp)./2.^(k-2);
Y(2,5) = x(5) -  (x(5)*perturbation_amp)./2.^(k-2);
Y(3,5) = x(5);

Y(1,6) = x(6) +  (x(6)*perturbation_amp)./2.^(k-2);
Y(2,6) = x(6) -  (x(6)*perturbation_amp)./2.^(k-2);
Y(3,6) = x(6);

Y(1,7) = x(7) +  (x(7)*perturbation_amp)./2.^(k-2);
Y(2,7) = x(7) -  (x(7)*perturbation_amp)./2.^(k-2);
Y(3,7) = x(7);

Y(1,8) = x(8) +  (x(8)*perturbation_amp)./2.^(k-2);
Y(2,8) = x(8) -  (x(8)*perturbation_amp)./2.^(k-2);
Y(3,8) = x(8);

Y(1,9) = x(9) +  (x(9)*perturbation_amp)./2.^(k-2);
Y(2,9) = x(9) -  (x(9)*0.2)./2.^(k-2);
Y(3,9) = x(9);

Y(1,10) = x(10) +  (x(10)*perturbation_amp)./2.^(k-2);
Y(2,10) = x(10) -  (x(10)*perturbation_amp)./2.^(k-2);
Y(3,10) = x(10);

Y(1,11) = x(11) +  (x(11)*perturbation_amp)./2.^(k-2);
Y(2,11) = x(11) -  (x(11)*perturbation_amp)./2.^(k-2);
Y(3,11) = x(11);

end

function [var1,var2,var3,var4,var5,var6,var7,var8,var9,var10,var11] = parameter_Assigning(x,M,index_1,index_2)

x(index_1) = M(index_2,index_1);
var1 = x(1);
var2 = x(2);
var3 = x(3);
var4 = x(4);
var5 = x(5);
var6 = x(6);
var7 = x(7);
var8 = x(8);
var9 = x(9);
var10 = x(10);
var11 = x(11);
end

function M = randomizer(x,n,m)
rng shuffle
M = zeros(m,n*length(x));
for j = 1:m
    temp = [];
    for i = 1:n
        temp = [temp x(randperm(length(x)))];
    end
    M(j,:) = temp;
end
end