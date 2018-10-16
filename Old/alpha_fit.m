clc
close all
clc
% a*(a*tanh(b*(x-c))+d)

f = 0:0.01:3.5;

a = 0.5107;
b = -1.377;
c = 1.968;
d = 1.467;

alpha = a*(a*tanh(b*(f-c))+d);

plot(f,alpha)