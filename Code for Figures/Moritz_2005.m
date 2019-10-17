
close all
clear all
clc

force_vec = [2 5 15 30 50 70 85 95];

CoV_vec = [6/7*0.5+4.5 ... 
    2/7*0.5+2.5 ...
    5/7*0.5+1 ...
    1.5 ...
    3.5/7*0.5 + 1.5 ...
    6/7*0.5+1 ...
    3/7*0.5+1 ...
    2/7*0.5+1];

figure(1)
plot(force_vec,CoV_vec)
hold on 
plot(force_vec,CoV_vec,'o')

%%
U = [0:0.01:1]*100;
U_th = 0.01*100;
CV_ISI = zeros(1,length(U));
for i = 1:length(U)
    CV_ISI(i) = 10+20*exp(-(U(i)-U_th)/2.5);
end

plot(U,CV_ISI)