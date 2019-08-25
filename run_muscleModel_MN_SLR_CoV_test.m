%==========================================================================
% run_muscleModel.m
% Author: Akira Nagamori
% Last update: 8/21/19
% Descriptions:
%   Run muscle model simulation
%==========================================================================
close all
clear all
clc


%%
data_folder = '/Volumes/DATA2/New_Model/SLR/CoV_ISI_test';
code_folder = '/Users/akira/Documents/Github/Twitch-Based-Muscle-Model';
model_parameter_folder =  '/Users/akira/Documents/Github/Twitch-Based-Muscle-Model/Model Parameters/Model_7';

%%
cd(model_parameter_folder)
load('modelParameter')
load('parameterMN')
cd(code_folder)
%% MU simulation parameters
modelParameter.CV_MU = 0.1;

%% Recruitment Type
modelParameter.recruitment = 3; % 1: Loeb's formulation, 2: Fuglevand's formulation


%% Simlulation parameters

amp_vec = [0.01:0.01:1];

for j = 31:length(amp_vec)
    j
    Fs = 30000;
    time = 0:1/Fs:15;

    %% 
    controlOpt = 1;
    % 1: feedfoward input
    % 2: feedback 
    
    %%
    SLRParameter.gamma_dynamic = 20;
    SLRParameter.gamma_static = 20;
    SLRParameter.Ia_delay = 20*Fs/1000;
    SLRParameter.Ia_gain = 10000;
    
    SLRParameter.G1 = 60; % conversion factor (Hz)
    SLRParameter.G2 = 4; % conversion factor (N)
    % transfer function describing GTO dynamics
    s = tf('s');
    H = (1.7*s^2+2.58*s+0.4)/(s^2+2.2*s+0.4);
    Hd = c2d(H,1/Fs);
    [num,den] = tfdata(Hd);
    SLRParameter.num_GTO = cell2mat(num);
    SLRParameter.den_GTO = cell2mat(den);
    SLRParameter.Ib_gain = 10000;
    SLRParameter.Ib_delay = 40*Fs/1000;
    
    delta = 0.0015;
    tau1 = 0.14;
    tau3 = 0.003;
    tau4 = 0.09;
    H_RI = (1+tau1*s)*exp(-delta*s)/((1+tau3*s)*(1+tau4*s));
    Hd_RI = c2d(H_RI,1/Fs);
    [num_RI_temp,den_RI_temp] = tfdata(Hd_RI);
    SLRParameter.num_RI = cell2mat(num_RI_temp);
    SLRParameter.den_RI = cell2mat(den_RI_temp);

    SLRParameter.RI_gain = 10;
    SLRParameter.RI_delay = 5*Fs/1000;
    SLRParameter.C_delay = 200*Fs/1000;
    SLRParameter.K_C = 0.0002;
    
    SLRParameter.noise_amp_Ia = 0;
    SLRParameter.noise_amp_Ib = 0;
    SLRParameter.noise_amp_RI = 0;
    SLRParameter.noise_amp_C = 0;
    SLRParameter.noise_amp_ID = 10000;
    SLRParameter.noise_amp_CD = 1;
    
    amp = amp_vec(j);
    input = [zeros(1,1*Fs) amp/2*[0:1/Fs:2] amp*ones(1,length(time)-1*Fs-length(amp*[0:1/Fs:2]))];
    %input_e =0.01*sin(2*pi*10*time);
    %input(5*Fs+1:end) = input(5*Fs+1:end) + input_e(5*Fs+1:end);
    %%
    
    for i = 1
        i
        tic
        output = spikeDrivenMuscleModel_MN_SLR(Fs,time,input,modelParameter,parameterMN,SLRParameter,controlOpt,1);
        toc
        cd(data_folder)
        save(['Data_' num2str(j) '_' num2str(i)],'output','-v7.3')
        cd(code_folder)
%         clear output

    end
    
end

