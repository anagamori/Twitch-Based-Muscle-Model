%--------------------------------------------------------------------------
% CT_HRT_fit_Stephens.m
% Author: Akira Nagamori
% Last update: 7/19/18
% Code descriptions
%   Find an equation that describes the length-dependence of contraction 
%   time (CT) and half-relaxation time (HRT) based on the data presented 
%   in Fig 4 A & B presented in Fig 9 C in Stephens et al. (1975)
%--------------------------------------------------------------------------

close all
clear all
clc

CT = [73 74 75.5 76.5 78.5 82 86.5 91 96.5 100 105 107.5 109.5 109]/100;
Lce_raw = -18:2:8; 
fiber_length = 3*10;
Lce = (fiber_length+Lce_raw)/fiber_length;

a1 = 0.416;
b1 = 1.235;
c1 = 0.45;
a2 = 1.015;
b2 = -9.718;
c2 = 17.24;
       
Lce_fit = linspace(Lce(1),Lce(end),length(CT)*5);
CT_fit =  a1*exp(-((Lce_fit-b1)/c1).^2) + a2*exp(-((Lce_fit-b2)/c2).^2);

figure(1)
plot(Lce,CT,'o','LineWidth',1)
hold on 
plot(Lce_fit,CT_fit,'LineWidth',1)

%%
HRT = [68 68 69 71 73 77.5 82 87.5 94 102 111 118 124 126]/100;
a1 = 1.114;
b1 = 1.446;
c1 = 0.7234;
a2 = 0.09588;
b2 = 1.192;
c2 = 0.194;
a3 = 0.8869;
b3 = -0.4906;
c3 = 1.271;
HRT_fit =  a1*exp(-((Lce_fit-b1)/c1).^2) + a2*exp(-((Lce_fit-b2)/c2).^2) + a3*exp(-((Lce_fit-b3)/c3).^2);

figure(2)
plot(Lce,HRT,'o','LineWidth',1)
hold on 
plot(Lce_fit,HRT_fit,'LineWidth',1)
       
