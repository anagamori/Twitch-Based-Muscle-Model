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
K_1 = 1.4;
K_2 = 0.13;
T = 350;
%--------------------------------------------------------------------------
% initialization
y_0 = 0; % detached state
y_1 = 0; % attached state

%--------------------------------------------------------------------------
% simulation parameters
Fs = 1000;
time = 0:1/Fs:5; %simulation time

stim = zeros(1,length(time));


f = f_0*(c^2/(c^2+c*K_2+K_1*K_2));
dy_1 = f*y_0 - g*y_1;