%close all
%clear all
%clc

testUnit = 184;
U_vec = 0:0.001:1;

U_th = U_th_new(testUnit); %0.01;
f_half = MDR(testUnit); %8;
f_2 = f_half*4;
f_onehalf = f_half*3;
alpha = 30;

a = (2*f_half+alpha*f_half)/(alpha*(1-U_th));
x = (a-f_half)/a;

tau = 0.01;

for i = 1:length(U_vec)
    U = U_vec(i);
    if U <= x
        DR(i) = f_half + alpha*a*(U-U_th);
    else
        DR(i) = f_2-a*(1-U);
    end
    if U < tau
        lamda(i) = 0;
    else
        lamda(i) = 19+8*U-(21+116*exp(-U/0.2))*tau...
            -exp((tau-U)/(0.16*tau+0.04))*(9.9+8*U+(-14.7-116*exp(-U/0.2))*tau);    
    end
end

DR(DR<f_half) = 0;

figure(1)
plot(U_vec,DR_mat(testUnit,:)./DR_mat(testUnit ,end))
hold on
plot(U_vec,DR./DR(end))
hold on
plot(U_vec,lamda./lamda(end))
