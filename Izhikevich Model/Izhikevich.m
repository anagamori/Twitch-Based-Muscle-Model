%==========================================================================
% Izhikevich_2003.m
% Author: Akira Nagamori
% Last update: 6/18/19
% Reference:
%   Izhikevich 2003 
%==========================================================================


function [v_vec,binary] = Izhikevich(time,input,parameter,Fs)

%%
I = input;
%% Model parameters
a = parameter.a; %0.005; %[0.02*ones(Ne,1); 0.02+0.08*ri];
b = parameter.b; %0.2; %[0.2*ones(Ne,1); 0.25-0.05*ri];
c = parameter.c; %-65; %[-65+15*re.^2; -65*ones(Ni,1)];
d = parameter.d; %8; %[8-6*re.^2; 2*ones(Ni,1)];

alpha = parameter.alpha; % 0.04
beta = parameter.beta; % 5
gamma = parameter.gamma; % 140

%% Initial Conditions
v = parameter.v; %-65; %*ones(Ne+Ni,1); % Initial values of v
u = b.*v; % Initial values of u

%%
binary = zeros(1,length(time));
v_vec = zeros(1,length(time));
for t = 1:length(time) 

    if v >= 30
    binary(t) = 1;
    v = c;
    u = u +d;
    end

v_dot = (alpha*v.^2+beta*v+gamma-u+I(t)); 
v = v_dot*1000/Fs + v;

u_dot = a.*(b.*v-u);
u = u_dot*1000/Fs + u; 

v_vec(t) = v;
end

end
