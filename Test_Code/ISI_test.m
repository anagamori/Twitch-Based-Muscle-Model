close all
clear all
clc

rng shuffle 
cv = 0.1;
DR = 40;
mu = 1/DR;
Z = randn(1,1000);
hist_ISI = mu + mu*cv*Z;

hist_ISI_2 = hist_ISI;
hist_ISI_2(hist_ISI_2<0.02) = 0.02;

figure(1)
histogram(hist_ISI)
hold on 
histogram(hist_ISI_2)

std(hist_ISI)/mean(hist_ISI)

std(hist_ISI_2)/mean(hist_ISI_2)