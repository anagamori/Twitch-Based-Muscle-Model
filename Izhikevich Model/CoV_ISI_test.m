%==========================================================================
% CoV_ISI_test.m
% Author: Akira Nagamori
% Last update: 8/21/19
%==========================================================================

close all
clear all
clc

code_folder = '/Users/akiranagamori/Documents/Github/Twitch-Based-Muscle-Model/Izhikevich Model';
model_parameter_folder =  '/Users/akiranagamori/Documents/Github/Twitch-Based-Muscle-Model/Model Parameters/Model_7';

%%
Fs = 10000;
time = 0:1/Fs:5;

U_amp_vec = 0.1; %0.01:0.01:1;

FR = zeros(1,length(U_amp_vec));
ISI_CoV = zeros(1,length(U_amp_vec));
r_vec = zeros(1,length(U_amp_vec));
%%
n = 1;
cd(model_parameter_folder)
load('modelParameter')
cd(code_folder)
FR_half = modelParameter.FR_half;
MDR = modelParameter.FR_half/2;
PDR = modelParameter.FR_half*2;
U_th = modelParameter.U_th;
testUnit = n;

U_th_target = U_th(testUnit);
MDR_target = MDR(testUnit);
PDR_target = PDR(testUnit);
FR_half_target = FR_half(testUnit);

cd([model_parameter_folder '/MN'])
load(['MN_' num2str(n)])
cd(code_folder)
    
%%
for i = 1:length(U_amp_vec)
    i
    U_amp = U_amp_vec(i);
    input = [zeros(1,1*Fs) U_amp*[0:1/Fs:2] U_amp*ones(1,length(time)-1*Fs-length(U_amp*[0:1/Fs:1]))];
    
    I_initial = (parameter.I_max-parameter.I_th)./(1-U_th(testUnit)).*(0-U_th(testUnit)) + parameter.I_th;
    v = (-5+0.2-sqrt((5-0.2)^2-4*0.04*(140+I_initial)))/(2*0.04);
    u = 0.2*v; 
    
    noise_1 = zeros(1,1);
    noise_2 = zeros(1,1);
     %length(time));
    v_vec = zeros(1,length(time));
    spike_train = zeros(1,length(time));
    I_vec = zeros(1,length(time));
    
    for t = 1:length(time)
        [noise_1] = noise(noise_1,Fs,10000);
        [noise_2] = noise(noise_2,Fs,1);
        U = input(t); %*(input(t));
        I = (parameter.I_max-parameter.I_th)/(1-U_th(testUnit))*(U-U_th(testUnit)) + parameter.I_th;
        I = I + I*(noise_1) + I*(noise_2);
        spike_vec = zeros(1,1);
        [u,v,spike_vec] = motoneuron(I,u,v,spike_vec,parameter.a,Fs);
        spike_train(t) = spike_vec;
        v_vec(t) = v;
        I_vec(t) = I; 
        
    end

    spike_time = find(spike_train(3*Fs+1:end));
    ISI = diff(spike_time)/(Fs/1000);
    FR(i) = mean(1./ISI*1000);
    FR_sd = std(1./ISI*1000);
    ISI_CoV(i) = FR_sd/FR(i)*100
    
    temp1 = ISI(1:end-1);
    temp2 = ISI(2:end);
    [r,p] = corrcoef(temp1,temp2);
    r_vec(i) = r(1);
        
end

%%
% offset = mean(ISI_CoV(end-10:end));
% [peak,loc] = max(ISI_CoV);
% 
% dForce = (U_amp_vec-U_amp_vec(loc))*100;
% dForce(dForce<=0) = 0;
% CoV_ISI_ex = offset+(peak-offset)*exp(-dForce/2.5);
% 
% figure(1)
% plot(U_amp_vec*100,ISI_CoV,'LineWidth',1)
% hold on 
% plot(U_amp_vec*100,CoV_ISI_ex,'LineWidth',1)
% plot([U_th_target*100 U_th_target*100],[0 peak+5],'--k','LineWidth',1)
% legend('New model','Mortiz et al. 2005','Recruitment Threshold')
% xlabel('Actovaton (%Maximum)','FontSize',14)
% ylabel('CoV ISI (%)','FontSize',14)
%%
% figure(2)
% plot(U_amp_vec*100,FR,'LineWidth',1)
% hold on
% plot([U_th_target*100 U_th_target*100],[0 max(FR)+5],'--k','LineWidth',1)
% plot([0 100],[MDR_target MDR_target],'--r','LineWidth',1)
% plot([0 100],[PDR_target PDR_target],'--r','LineWidth',1)
% xlabel('Actovaton (%Maximum)','FontSize',14)
% ylabel('Discharge Rate (Hz)','FontSize',14)
%%
function [u,v,spike_vec] = motoneuron(I,u,v,spike_vec,a,Fs)
b_MN = 0.2;
c_MN = -65;
d_MN = 6;

alpha_MN = 0.04;
beta_MN  = 5;
gamma_MN = 140;

V_thres = 30+normrnd(0,1);
index = v>=V_thres;
spike_vec(index) = 1;
v(index) = c_MN+normrnd(0,1);
u(index) = u(index) + d_MN;

v_dot = (alpha_MN*v.^2+beta_MN.*v+gamma_MN-u+I);
v = v_dot*1000/Fs + v;

u_dot = a.*(b_MN.*v-u);
u = u_dot*1000/Fs + u;

end

function [x] = noise(x,Fs,amp)
    D = amp;
    tau = 0.005; 
    chi = normrnd(0,1,[1,size(x,2)]);
    x_dot = -x/tau + sqrt(D)*chi;
    x = x_dot*1/Fs + x;
end