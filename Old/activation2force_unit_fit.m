close all
clc

testintUnit = 1;
x = 8:8:48;
%y = [1, 1, 0.7, 0.64,0.6,0.56,0.5,0.45,0.41,0.37,0.34,0.31];
% y = [1, 1, 1, 0.672, 0.63,0.61,0.547,0.465,0.41,0.35,0.305,0.265];
y = [1 0.68 0.62 0.47 0.35 0.267];
f = x./FR_half(testintUnit);
f_test = 0:0.01:3;

fit1 = fit(f',y','poly1');
figure(1)
plot(fit1,f,y)

%%
a_1 = 0.69;
b_1 = 5;
alpha_1_temp = 0.28*(1-exp(-(f_test/a_1).^b_1))+1; 
alpha_2_temp = 2.031*exp(-0.4015*f_test) - 52.85*exp(-6.491*f_test);
%alpha_2_temp = 0.28*(1-exp(-(f_test/a_2).^b_2))+1; 
figure(2)
plot(f,y,'o')
% hold on
% plot(f_test(1:100),alpha_1_temp(1:100))
% hold on
% plot(f_test(100:end),alpha_2_temp(100:end))

