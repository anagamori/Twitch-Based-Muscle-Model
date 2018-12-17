%==========================================================================
% test_y_z.m
% Author: Akira Nagamori
% Last update: 12/17/18
% Descriptions:
%   Find a value of y (0-1) to get a desired steady-state value of
%   cross-bridge activation, z. 
%==========================================================================

close all
clear all
clc

% The dynamics of cross-bridge activation (i.e., fraction of attached cross bridge) 
% is defined as dz/dt = f_app*(1-z) - g_app*z (Eq.1) following the original formulation by Huxley  (Brenner, 1988)
% Then the steady-state value of z, z_ss, is given by z_ss = f_app/(f_app+g_app) (Eq.2);

% The apparant rate costant for cross-bridge detachment is constant across
% various levels of calcium (Sweeney & Stull, 1990)
g_app = 10; 

% Motoneuron firing rate 
FR = 0:0.01:2; %[f_0.5]
% Activation-frequency relationship from Song et al. (2008)
Af_Song = 1-exp(-(FR./(0.56*2.1)).^2.1);

%--------------------------------------------------------------------------
% find a value of the desired f_app that satisfies the activation-frequency
% relationship prposed by Song et al. (2008)
f_app_desired = zeros(1,length(FR));
for i = 1:length(FR)
    % Use Eq.2 to find a value of f_app given g_app and z_ss
    z_ss = sqrt(Af_Song(i));
    f_app_desired(i) = z_ss*g_app/(1-z_ss);
end

f_app_desired_norm = f_app_desired/f_app_desired(end);
%
figure(1)
plot(FR,f_app_desired,'LineWidth',2)
xlabel('FR','FontSize',14)
ylabel('Desired f_{app}','FontSize',14)

%--------------------------------------------------------------------------
% Fit a exponential curve to f_app to FR relationship
%--------------------------------------------------------------------------
% the fraction of available cross-bridges, y, ranges from 0 to 1, where y
% reaches the value of 1 when FR = 2
y = FR/2; 
f1 = fit(y',f_app_desired_norm','exp1');
figure(2)
plot(f1,y,f_app_desired_norm)
xlabel('FR','FontSize',14)
ylabel('Normalized Desired f_{app}','FontSize',14)



