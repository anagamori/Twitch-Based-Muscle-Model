
close all
clear all
clc

tau = 0.01;
phi_vec = 0:0.001:1;

for i = 1:length(phi_vec)
    phi = phi_vec(i);
    lamda(i) = 19+8*phi-(21+116*exp(-phi/0.2))*tau...
        -exp((tau-phi)/(0.16*tau+0.04))*(9.9+8*phi+(-14.7-116*exp(-phi/0.2))*tau);
end

figure(1)
plot(phi_vec,lamda)