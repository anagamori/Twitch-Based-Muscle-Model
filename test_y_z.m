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
% is defined as dz/dt = f_app*(1-z) - g_app*z following the original formulation by Huxley  (Brenner, 1988)
% Then the steady-state value of z, z_ss, is given by 
z_ss = f_app/(f_app+g_app);