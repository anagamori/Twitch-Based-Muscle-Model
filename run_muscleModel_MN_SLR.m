%==========================================================================
% run_muscleModel.m
% Author: Akira Nagamori
% Last update: 3/5/19
% Descriptions:
%   Run muscle model simulation
%==========================================================================
close all
clear all
clc


%%
data_folder = '/Volumes/DATA2/New_Model/SLR/GD_40_GS_40_Ia_2000_DL_30';
code_folder = '/Users/akira/Documents/Github/Twitch-Based-Muscle-Model';
model_parameter_folder =  '/Users/akira/Documents/Github/Twitch-Based-Muscle-Model/Model Parameters/Model_5_constantT2T';

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

amp_vec = [0.05 0.1:0.1:1];
trial_vec = [7 10];
for j = 1
    j
    if j < 2
        Fs = 10000;
        time = 0:1/Fs:10;
    elseif j >= 2 && j < 4
        Fs = 15000;
        time = 0:1/Fs:15;
    elseif j >= 4 && j < 7
        Fs = 20000;
        time = 0:1/Fs:15;
    elseif j >= 7
        Fs = 25000;
        time = 0:1/Fs:15;
    end
    %% 
    controlOpt = 1;
    % 1: feedfoward input
    % 2: feedback 
    %%
    SLRParameter.gamma_dynamic = 20;
    SLRParameter.gamma_static = 20;
    SLRParameter.Ia_delay = 20*Fs/1000;
    SLRParameter.Ia_gain = 10000;
    SLRParameter.Ib_gain = 10000;
    SLRParameter.Ib_delay = 40*Fs/1000;
    SLRParameter.RI_gain = 10;
    SLRParameter.RI_delay = 5*Fs/1000;
    SLRParameter.C_delay = 200*Fs/1000;
    SLRParameter.K_C = 0.0001;
    
    amp = amp_vec(j+1);
    input = [zeros(1,1*Fs) amp/2*[0:1/Fs:2] amp*ones(1,length(time)-1*Fs-length(amp*[0:1/Fs:2]))];
    %input_e =0.01*sin(2*pi*10*time);
    %input(5*Fs+1:end) = input(5*Fs+1:end) + input_e(5*Fs+1:end);
    %%
    
    for i = 1
        i
        tic
        output = spikeDrivenMuscleModel_MN_SLR(Fs,time,input,modelParameter,parameterMN,SLRParameter,controlOpt,1);
        toc
%         cd(data_folder)
%         save(['Data_' num2str(j) '_' num2str(i)],'output','-v7.3')
%         cd(code_folder)
%         clear output

    end
    
end
%%
temp_1 = output.ForceTendon(5*Fs+1:end);
[pxx,f] = pwelch(temp-mean(temp),gausswin(5*Fs),5*Fs*0.9,0:0.1:100,Fs,'power');
figure(2)
plot(f,pxx,'LineWidth',2)
xlim([0 30])
xlabel('Frequency (Hz)','FontSize',14)
ylabel('Power (N^2)','FontSize',14)

Force_CoV = std(temp)/mean(temp)*100

%%
ISI_vec = [];
ISI_sd_vec = [];
ISI_CoV_vec = [];
r_vec = [];
p_vec = [];
for n = 1:300
    index = 0;
    spike_time = find(output.spike_train(n,5*Fs+1:end));
    if ~isempty(spike_time)       
        ISI = diff(spike_time)/(Fs/1000);
        if length(ISI) >= 10
            
            FR_vec = nanmean(1./ISI*1000);
            FR_sd = nanstd(1./ISI*1000);
            ISI_CoV_vec = [ISI_CoV_vec FR_sd/FR_vec*100];
            ISI_sd_vec = [ISI_sd_vec nanstd(ISI)];
            ISI_vec = [ISI_vec nanmean(ISI)];
            temp1 = ISI(1:end-1);
            temp2 = ISI(2:end);
            [r,p] = corrcoef(temp1,temp2);
            r_vec = [r_vec r(1,2)];
            p_vec = [p_vec p(1,2)];
        end
    end
end
figure(5)
histogram(ISI_CoV_vec,[0:5:100])

%%
figure()
histogram(r_vec(p_vec<0.05),-0.5:0.05:0.5)
hold on 
histogram(r_vec(p_vec>=0.05),-0.5:0.05:0.5)
%%
spike_time = find(output.spike_train(6,:));
FR_vec = zeros(1,length(time));
for i = 2:length(spike_time)
    ISI_temp = (spike_time(i) - spike_time(i-1))/(Fs/1000);
    FR = 1./ISI_temp*1000;
    FR_vec(spike_time(i-1):spike_time(i)) = FR;
    if i == 2
        FR_vec(spike_time(i-1):spike_time(i)) = FR;
    elseif i == length(spike_time)
        FR_vec(spike_time(i):end) = FR;
    end
end
figure()
plot(time,FR_vec)
xlabel('Time (s)','FontSize',14)
ylabel('Discharge Rate (Hz)','FontSize',14)
%% 
temp1 = zeros(1,length(time));
temp1(find(output.spike_train(1,:))) = 1;
temp1 = temp1(5*Fs+1:end);
temp2 = zeros(1,length(time));
temp2(find(output.spike_train(10,:))) = 1;
temp2 = temp2(5*Fs+1:end);
[C,f] = mscohere(temp1-mean(temp1),temp2-mean(temp2),gausswin(2*Fs),2*Fs*0.9,0:0.5:100,Fs);
figure(4)
plot(f,C)
xlabel('Frequency (Hz)','FontSize',14)
ylabel('Cohrence','FontSize',14)