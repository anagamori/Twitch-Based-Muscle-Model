%==========================================================================
% Todorov_2005.m
% Author: Akira Nagamori
% Last update: 3/1/20
% Model descriptions:
%   input: Fs = sampling frequency
%             time = time vector
%             synaptic_drive = synaptic drive to a motor unit population
%             modelParameter = a structure that contains all model
%             parameters
%             fitOpt = figure display option (0: no figure output, 1: display
%             figures)
%   output: a data structure ('output') that contains spike trains, motor
%   unit forces and tendon force (see l.242)
%==========================================================================

Fs = 100; % sampling frequency
time = 0:1/Fs:15; % time vector

amp_vec = 0.1:0.1:1; % 

% initialization
mean_f = zeros(1,length(amp_vec));
std_f = zeros(1,length(amp_vec));

for i = 1:length(amp_vec)
    rng shuffle
    
    u = zeros(1,length(time));
    u(3*Fs+1:end) = amp_vec(i);
    r = normrnd(0,1,[1,length(time)]);
    omega_c = 0.103;
    
    tau_1 = 0.03;
    tau_2 = 0.04;
    
    mean_f_temp = zeros(1,20);
    std_f_temp = zeros(1,20);
    
    for j = 1:10
        f = 0;
        f_dot = 0;
        
        u_vec = zeros(1,length(time));
        f_vec = zeros(1,length(time));
        
        for t = 1:length(time)
            u_t = u(t)*omega_c*r(t)+u(t);
            f_ddot = (u_t - f - (tau_1+tau_2)*f_dot)/(tau_1*tau_2);
            f_dot = f_ddot/Fs + f_dot;
            f = f_dot/Fs + f;
            
            u_vec(t) = u_t;
            f_vec(t) = f;
        end
        
        mean_f_temp(j) = mean(f_vec(3*Fs+1:end));
        std_f_temp(j) = std(f_vec(3*Fs+1:end));
        
    end
    mean_f(i) = mean(mean_f_temp);
    std_f(i) = mean(std_f_temp);
    
end

figure(2)
plot(mean_f,std_f)

f = polyfit(mean_f,std_f,1)