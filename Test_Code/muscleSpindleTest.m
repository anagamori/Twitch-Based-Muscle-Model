% close all
% clear all
% clc

Fs = 10000; % sampling frequency of simulation
t = 0:1/Fs:5; % time vector

% simulated length change of muscle
% First column of Fig 3 in Mileusnic et al., (2006)
lengthChange = [0.95.*ones(1,1*Fs) 0.11*t(1:1.1818*Fs)+0.95 1.08*ones(1,length(t)-1*Fs-1.1818*Fs)];
% smooth out the length change to avoid sharp transitions, which might
% affect differentiation 
Lce = output.Lce;  %smooth(lengthChange,100);
%Lce = Lce'; % muscle length
Vce = output.Vce; %(2:end); %[0 diff(Lce).*Fs]; % muscle velocity
Ace = output.Ace; %(2:end); %[0 diff(Vce).*Fs]; % muscle acceleration


% Fusimotor drive
gamma_dynamic = 20; % gamma dynamic
gamma_static = 20; % gamma static

% Run a simulation 
[outputIa,outputII,AP_bag1,AP_primary_bag2,AP_primary_chain] = muscleSpindle(Lce,Vce,Ace,Fs,gamma_dynamic,gamma_static);

