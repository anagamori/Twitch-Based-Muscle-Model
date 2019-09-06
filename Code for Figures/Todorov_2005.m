close all
clear all
clc

Fs = 100;
time = 0:1/Fs:10;

amp_vec = 0.1:0.1:1;

mean_f = zeros(1,length(amp_vec));
std_f = zeros(1,length(amp_vec));

for i = 1:length(amp_vec)
rng shuffle

u = zeros(1,length(time));
u(2*Fs+1:end) = amp_vec(i);
r = normrnd(0,1,[1,length(time)]);
omega_c = 0.4;

tau_1 = 0.04;
tau_2 = 0.04;

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

mean_f(i) = mean(f_vec(3*Fs+1:end));
std_f(i) = std(f_vec(3*Fs+1:end));

figure(1)
subplot(2,1,1)
plot(time,u)
hold on 
plot(time,u_vec)
subplot(2,1,2)
plot(time,f_vec)
hold on

end

figure(2)
plot(mean_f,std_f)