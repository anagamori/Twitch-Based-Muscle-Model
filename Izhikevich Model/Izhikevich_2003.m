%==========================================================================
% Izhikevich_2003.m
% Author: Akira Nagamori
% Last update: 6/18/19
% Reference:
%   Izhikevich 2003 
%==========================================================================

close all
clear all
clc

%%
Fs = 10000;
time = 0:1/Fs:5;
%%
I = zeros(1,length(time));
I(2*Fs+1:end) = 4;
%% Model parameters
a= 0.005; %[0.02*ones(Ne,1); 0.02+0.08*ri];
b= 0.2; %[0.2*ones(Ne,1); 0.25-0.05*ri];
c= -65; %[-65+15*re.^2; -65*ones(Ni,1)];
d= 8; %[8-6*re.^2; 2*ones(Ni,1)];

%% Initial Conditions
v=-65; %*ones(Ne+Ni,1); % Initial values of v
u=b.*v; % Initial values of u

%%
binary = zeros(1,length(time));
v_vec = zeros(1,length(time));
for t = 1:length(time) 

    if v >= 30
    binary(t) = 1;
    v = c;
    u = u +d;
    end

v_dot = (0.04*v.^2+5*v+140-u+I(t)); 
v = v_dot*1000/Fs + v;

u_dot = a.*(b.*v-u);
u = u_dot*1000/Fs + u; 

v_vec(t) = v;
end

figure(1)
subplot(2,1,1)
plot(time,I)
subplot(2,1,2)
plot(time,v_vec)

spike_time = find(binary(3*Fs+1:end));
ISI = diff(spike_time)/(Fs/1000);
mean_FR = mean(1./ISI*1000)

