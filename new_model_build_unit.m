%==========================================================================
% new_model_build_unit.m
% Author: Akira Nagamori
% Last update: 2/24/119
%==========================================================================

close all
clear all
clc

% Create a templatee using new_model_test.m
% Take its parameter and use new_model_parameterFit to fine-tune the
% parameters
new_model_parameterFit_ST
% Use the fitted model to then fit it to length-dependence of
% activation-frequency relationship using new_model_parameterFit_ST_length
new_model_parameterFit_ST_length
