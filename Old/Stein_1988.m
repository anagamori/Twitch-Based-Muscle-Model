%==========================================================================
% Stein_1988.m
% Author: Akira Nagamori
% Last update: 10/29/18
% Descriptions: find parameters, a and b, for non-linear scaling of force
% output, for each motor unit
%==========================================================================
close all
clear all
clc
%--------------------------------------------------------------------------
% constant parameters
f_0 = 3;
g_0 = 2.1;
g_1 = 13.5;
g = g_0;
K_1 = 1.4;
K_2 = 0.13;
T = 350;
am = 23.4;
b = 3.6;
p = b/(am);
q = 1/(am);
%--------------------------------------------------------------------------
% initialization
y_0 = 1; % detached state
y_1 = 0; % attached state
dF = 0;
F = 0;
trigger = 0;
%--------------------------------------------------------------------------
% simulation parameters
Fs = 1000;
time = 0:1/Fs:3; %simulation time

stim = zeros(1,length(time));
stim(1*Fs) = 1; 

for t = 1:length(time)
    if stim(t) == 1
        trigger = t;
    end
    
    if t-trigger < T
        g = g_0;
    else
        g = g_1;
    end
    
    %c = p*F + q*dF;
    c = 0.59;
    f = f_0*(c^2/(c^2+c*K_2+K_1*K_2));
    dy_1 = f*y_0*stim(t) - g*y_1;
    y_1 = dy_1*(1/Fs) + y_1;
    
    F_r = y_1*(f_0+g)/f_0;
    
    y_1_vec(t) = y_1;
    F_r_vec(t) = F_r;
    
end

figure(1)
subplot(2,1,1)
plot(time,y_1_vec)
subplot(2,1,2)
plot(time,F_r_vec)




