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
% run multiple times at each length to fit a function 


%%==========================================================================
% To do list
%%==========================================================================
% 1. Streamline model development
%       Make m-files /into functions (e.g., new_model_build_test.m into a function that can be called inside new_model_parameterFit_ST)
%       Code that can be for-looped and carried out without interuption to
%       generate a single MU
% 2. Find a function to fit length-dependence of parameters
%       Run fitting multiple times to see the pattern
% 3. Build a second template for slow-twitch MU
% 4. Test the fist-twitch template for length-dependence
% 5. Create a muscle with 300 MUs

      

