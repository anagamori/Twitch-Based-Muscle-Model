%==========================================================================
% test_frequency_response.m
% Author: Akira Nagamori
% Last update: 6/21/19
%==========================================================================

close all
clear all
clc

code_folder = '/Users/akiranagamori/Documents/Github/Twitch-Based-Muscle-Model/Izhikevich Model';
model_parameter_folder =  '/Users/akiranagamori/Documents/Github/Twitch-Based-Muscle-Model/Model Parameters/Model_4_Ur_50_constantT2T';

%%
Fs = 10000;
time = 0:1/Fs:20;

[b,a] = butter(4,100/(Fs/2),'low');
%%
for n = 300
    cd(model_parameter_folder )
    load('modelParameter')
    cd(code_folder)
    FR_half = modelParameter.FR_half;
    MDR = modelParameter.FR_half/2;
    PDR = modelParameter.FR_half*2;
    U_th = modelParameter.U_th_new;
    testUnit = n;
    
    U_th_target = U_th(testUnit);
    MDR_target = MDR(testUnit);
    PDR_target = PDR(testUnit);
    FR_half_target = FR_half(testUnit);
    
    cd([model_parameter_folder '/MN'])
    load(['MN_' num2str(n)])
    cd(code_folder)
    
    U_amp = 0.3;
    input = [zeros(1,1*Fs) U_amp/2*[0:1/Fs:2] U_amp*ones(1,length(time)-1*Fs-length(U_amp*[0:1/Fs:2]))];
    
    I_initial = (parameter.I_max-parameter.I_th)./(1-U_th(testUnit)).*(0-U_th(testUnit)) + parameter.I_th;
    v = (-5+0.2-sqrt((5-0.2)^2-4*0.04*(140+I_initial)))/(2*0.04);
    u = 0.2*v;
    
    f_vec = 1:1:30;
    FR_p2p = zeros(1,length(f_vec));
    for f = 1:length(f_vec)
        
        x_noise = 0;
        v_vec = zeros(1,length(time));
        spike_train = zeros(1,length(time));
        I_vec = zeros(1,length(time));
        x_noise_vec = zeros(1,length(time));
        
        time_sin = 0:1/Fs:15;
        sin_test = 0.05*sin(2*pi*f_vec(f)*time_sin);
        input(5*Fs+1:end) = input(5*Fs+1:end)+sin_test;
        
        for t = 1:length(time)
            [x_noise] = noise(x_noise,Fs);
            U = input(t); %+x_noise*input(t);
            I = (parameter.I_max-parameter.I_th)/(1-U_th(testUnit))*(U-U_th(testUnit)) + parameter.I_th;
            spike_vec = zeros(1,1);
            [u,v,spike_vec] = motoneuron(I,u,v,spike_vec,parameter.a,Fs);
            spike_train(t) = spike_vec;
            v_vec(t) = v;
            I_vec(t) = I;
            x_noise_vec(t) = x_noise;
        end
        
        spike_time = find(spike_train(5*Fs+1:end));
        ISI = diff(spike_time)/(Fs/1000);
        FR = 1./ISI*1000;
        FR_p2p(f) = max(FR) - min(FR);
        
        
    end
    
    FR_dB =10*log10(FR_p2p./0.05);
figure(1)
% ax1 = subplot(2,1,1);
plot(f_vec,FR_p2p,'LineWidth',2)
% ax2 = subplot(2,1,2);
% plot(time,v_vec,'LineWidth',2)
% linkaxes([ax1 ax2],'x')

end



function [u,v,spike_vec] = motoneuron(I,u,v,spike_vec,a,Fs)
b_MN = 0.2;
c_MN = -65;
d_MN = 6;

alpha_MN = 0.04;
beta_MN  = 5;
gamma_MN = 140;

index = v>=30;
spike_vec(index) = 1;
v(index) = c_MN;
u(index) = u(index) + d_MN;

v_dot = (alpha_MN*v.^2+beta_MN.*v+gamma_MN-u+I);
v = v_dot*1000/Fs + v;

u_dot = a.*(b_MN.*v-u);
u = u_dot*1000/Fs + u;

end

function [x] = noise(x,Fs)
D = 100000;
tau = 0.01;
chi = normrnd(0,1,[1,1]);
x_dot = -x/tau + sqrt(D)*chi;
x = x_dot*1/Fs + x;
end