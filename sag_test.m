%==========================================================================
% sag_test.m
% Author: Akira Nagamori
% Last update: 4/2/19
%==========================================================================
close all
clc
clear all
%%
Fs = 2000;
time = 0:1/Fs:5;
amp = 1;

%% Loeb's model
%f_eff = [zeros(1,1*Fs) amp/2*[0:1/Fs:2] amp*ones(1,length(time)-1*Fs-length(amp*[0:1/Fs:2]))];
f_eff = [zeros(1,1*Fs) amp*ones(1,length(time)-1*Fs)];

a_f = 0.56;
n_f = 2.1;
S_i = 0.96;
Af = zeros(1,length(time));
S_vec = zeros(1,length(time));
for t = 1:length(time)
    [S_i] = sag_function(S_i,f_eff(t),Fs);
    S_vec(t) = S_i;
    Af(t) = 1-exp(-((S_i*f_eff(t))/(a_f*n_f))^n_f);
end

figure(1)
plot(time,Af,'k')
hold on
%% New model
code_folder = '/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model';
data_folder = '/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Model Parameters/Model_1';
n = 300;
cd(data_folder)
load('FR_half')
load('pool_parameter_matrix')
cd(code_folder)

FR_half_test = FR_half(n);
FR_test = FR_half_test*amp;
spike = zeros(1,length(time));
% Generate spike train
temp =spikeTrainGenerator(0:1/Fs:4,Fs,FR_test);
spike(1*Fs+1:end-1) = temp(1:end-1);

parameter = parameter_Matrix(n,:);
Lce = 1;

S = parameter(1); %7;
C = parameter(2); %1.025;
k_1 = parameter(3); %14.625;
k_2 = parameter(4); %4.9375;
k_3 = parameter(5)*Lce + parameter(6); %17.41*Lce - 2.85;
k_4 = parameter(7)*Lce + parameter(8); %-7.67*Lce + 14.92;
tau_1 = parameter(9); %0.0051;
tau_2 = parameter(10); % 0.04;
N = parameter(11)*Lce + parameter(12); %-2.26*Lce + 4.20;
K = parameter(13)*Lce + parameter(14); %-0.044*Lce + 0.080;
alpha = parameter(15); %4.475;

c = 0; % free calcium concentration
cf = 0; % concentraction of calcium bound to troponin
A = 0; % muscle activation

S_new = 0;
k_2_new = 0;
c_vec = zeros(1,length(time));
cf_vec = zeros(1,length(time));
A_tilda_vec = zeros(1,length(time));
A_vec = zeros(1,length(time));

R_temp = exp(-time/tau_1);
R = zeros(1,length(time));
for t = 1:length(time)
    %% Stage 1
    spike_temp = zeros(1,length(time));
    % Calcium diffusion to sarcoplasm
    if spike(t) == 1
        spike_temp(t) = 1;
        temp = conv(spike_temp,R_temp*(1+2*A^alpha));
        R = R + temp(1:length(time));
    end
    %R = spike(t) + exp(-h/tau_1)*R; %*(1+3*A^alpha);
    
    
    if f_eff(t) < 0.1
        a_s = 20;
    else
        a_s = 0.96;
    end
    
    T_s = 0.015;
    S_new_dot = (a_s - S_new)./T_s;
    S_new = S_new_dot/Fs + S_new;
    
    if f_eff(t) < 0.1
        k_2_temp = k_2;
    else
        k_2_temp = k_2*1.96;
    end
    
    T_s = 0.015;
    k_2_new_dot = (k_2_temp - k_2_new)./T_s;
    k_2_new  = k_2_new_dot/Fs + k_2_new ;
    
    %%
    c_dot = k_1*(C-c-cf)*R(t) - k_2*c*(S-C+c+cf)-(k_3*c-k_4*cf)*(1-cf);
    cf_dot = (1-cf)*(k_3*c-k_4*cf);
    c = c_dot/Fs + c;
    cf = cf_dot/Fs + cf;
    
    %% Stage 2
    % Cooperativity and saturation
    if cf < 0
        cf_temp = 0;
    else
        cf_temp = cf*S_new;
    end
    A_tilda = cf_temp^N/(cf_temp^N+K^N);
    
    %% Stage 3
    % First-order dynamics to muscle activation, A
    %tau_2 = tau_2_0*(1-0.8*A_tilda)^2;
    A_dot = (A_tilda-A)/tau_2;
    A = A_dot/Fs + A;
    
    %% Store variables
    %x_vec(t) = x(t);
    %R_vec(t) = R;
    c_vec(t) = c;
    cf_vec(t) = cf;
    A_tilda_vec(t) = A_tilda;
    A_vec(t) = A;
    
    S_i_new_vec(t) = S_new;
    k_2_new_vec(t) = k_2_new;
end

figure(1)
plot(time,A_vec,'b')
hold on
% plot(time,A_vec.*S_vec,'r')
% hold on
%plot(time,A_vec_new)
%%
function [S] = sag_function(S,f_eff,Fs)

if f_eff < 0.1
    a_s = 1.76;
else
    a_s = 0.96;
end

T_s = 0.043;
S_dot = (a_s - S)./T_s;
S = S_dot/Fs + S;

end

function spikeTrain = spikeTrainGenerator(t,Fs,freq)

spikeTrain = zeros(1,length(t));
ISI = round(1/freq*Fs);
numSpikes = round(length(t)/ISI);
index = [1:numSpikes]*ISI;
index(index>length(t)) = [];
spikeTrain(index) = 1;
spikeTrain(1) = 1;

end

function [var1,var2,var3,var4,var5,var6,var7,var8,var9,var10,var11] = parameter_Assigning(x)
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