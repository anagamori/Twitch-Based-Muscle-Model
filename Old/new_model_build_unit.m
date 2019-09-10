%==========================================================================
% new_model_build_unit.m
% Author: Akira Nagamori
% Last update: 3/9/119
%==========================================================================

close all
clear all
clc

%% Create a template
% Create a templatee using new_model_test.m
new_model
%% Fine-tune parameters
% Take its parameter and use new_model_parameterFit to fine-tune the
% parameters
new_model_parameterFit_ST 
new_model_parameterFit_FT
%% Fit length-dependence 
% Use the fitted model to then fit it to length-dependence of
% activation-frequency relationship using new_model_parameterFit_ST_length
% run multiple times at each length to fit a function 
new_model_parameterFit_ST_length
new_model_parameterFit_FT_length

% analyze the realtionship between length and a value of each parameter
model_Af_length_analysis 

%% Run a new model 
run_new_model



      

