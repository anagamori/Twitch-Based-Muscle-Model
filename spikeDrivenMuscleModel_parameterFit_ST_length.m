%==========================================================================
% new_model_parameterFit_ST_length.m
% Author: Akira Nagamori
% Last update: 2/22/119
%==========================================================================

close all
clear all
clc

%% Folder name
code_folder = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model';
data_folder = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Model Parameters/Model_1/ST';
%% Simulation parameters
Fs = 1000; %sampling frequency
T = 1/Fs;
time = 0:1/Fs:5; %simulation time

%% Parameters to be searched
%Lce = 1.1;
first_MU = 30;
last_MU = 196;
parpool(12)
for trialN = first_MU:last_MU
    trialN
    cd(data_folder)
    load(['Data_' num2str(trialN)])
    cd(code_folder)
    
    FR_half_initial = Data{2,6};
    param_initial = Data{2,12};
    
    Data_temp = cell(1,40);
    parfor j = 1:40
        
        if j <= 10
            Lce = 0.8;
        elseif j > 10 && j <= 20
            Lce = 0.9;
        elseif j > 20 && j <= 30
            Lce = 1.1;
        elseif j > 30
            Lce = 1.2;
        end
%         if j == 1
%             Lce = 0.8;
%         elseif j == 2
%             Lce = 0.9;
%         elseif j == 3
%             Lce = 1.1;
%         elseif j == 4
%             Lce = 1.2;
%         end
        condition = j+10;
        FR_half = FR_half_initial;
        param = param_initial;
        for k = 1:6
            rng shuffle
            Param_matrix = annealing_curve(param,k);
            
            r = randperm(4);
            
            %% Loop through all parameters
            for n = 1:length(r)
                if r(n) == 1
                    index = 5;
                elseif r(n) == 2
                    index = 6;
                elseif r(n) == 3
                    index = 9;
                elseif r(n) == 4
                    index = 10;
                end
                %% Loop through all perturbations
                error_long = zeros(1,3);
                for l = 1:3
                    [S,C,k_1,k_2,k_3,k_4,tau_1,tau_2,N,K,alpha] = parameter_Assigning(param,Param_matrix,index,l);
                    param =  [S,C,k_1,k_2,k_3,k_4,tau_1,tau_2,N,K,alpha];
                    %% Run a twitch simulation and sweep simulation
                    [Data_temp_2] = spikeDrivenMuscleModel_testFunction(param,Lce,FR_half,'slow');
                    
                    error_long(l) = Data_temp_2{2,8};
                    [min_error,loc_min_error] = min(error_long);
                end
                [S,C,k_1,k_2,k_3,k_4,tau_1,tau_2,N,K,alpha] = parameter_Assigning(param,Param_matrix,index,loc_min_error);
                param =  [S,C,k_1,k_2,k_3,k_4,tau_1,tau_2,N,K,alpha];
            end
        end
        
        [Data_temp{j}] = spikeDrivenMuscleModel_testFunction(param,Lce,FR_half,'slow');
        
    end
    for trial = 1:40
        Data  = Data_temp{trial};
        cd(data_folder)
        save(['Data_' num2str(trialN) '_' num2str(trial)],'Data')
        cd(code_folder)
    end
    close all
end
%%
function spikeTrain = spikeTrainGenerator(t,Fs,freq)

spikeTrain = zeros(1,length(t));
ISI = round(1/freq*Fs);
numSpikes = round(length(t)/ISI);
index = [1:numSpikes]*ISI;
index(index>length(t)) = [];
spikeTrain(index) = 1;
spikeTrain(1) = 1;

end

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

function error = error_calculation(vec_1,vec_2,reference)
for i = 1:length(find(reference<2.5&reference>0.4))
    error(i) = abs(vec_1(i)-vec_2(i));
end
end