%--------------------------------------------------------------------------
% Tw2Tet_fit_Stephens.m
% Author: Akira Nagamori
% Last update: 7/19/18
% Code descriptions
%   Find an equation that describes the length-dependence of twitch-tetanus ratio 
%   based on the data presented in Fig 9 C in Stephens et al. (1975)
%--------------------------------------------------------------------------

close all
clear all
clc

Tw2Tet = [0.11 0.122 0.145 0.184 0.224 0.245 0.264 0.272 0.275 0.272 0.261 0.244 0.22 0.189];
Tw2Tet = Tw2Tet./max(Tw2Tet);
Lce_raw = -18:2:8; 
fiber_length = 3*10;
Lce = (fiber_length+Lce_raw)/fiber_length;

% Guassian 4th order
a1 = 0;
b1 = 0.9607;
c1 = 0.003597;
a2 = 0.1882;
b2 = 0.706;
c2 = 0.1679;
a3 = 0.5776;
b3 = 0.9928;
c3 = 0.3403;
a4 = 0.4059;
b4 = 0.9397;
c4 = 1.541;

Lce_fit = linspace(Lce(1),Lce(end),length(Tw2Tet)*5);
Tw2Tet_fit = a1*exp(-((Lce_fit-b1)/c1).^2) + a2*exp(-((Lce_fit-b2)/c2).^2) +  a3*exp(-((Lce_fit-b3)/c3).^2) + a4*exp(-((Lce_fit-b4)/c4).^2);

figure(1)
plot(Lce,Tw2Tet,'o','LineWidth',1)
hold on 
plot(Lce_fit,Tw2Tet_fit,'LineWidth',1)