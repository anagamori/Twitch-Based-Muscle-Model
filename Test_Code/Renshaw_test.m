close all
clear all
clc

Fs = 10000;
time = 0:1/Fs:5;
amp = 0.1;
U = [zeros(1,1*Fs) amp/2*[0:1/Fs:2] amp*ones(1,length(time)-1*Fs-length(amp*[0:1/Fs:2]))];

U_mat = zeros(100,length(time));

s = tf('s');
delta = 0.0015;
tau1 = 0.14;
tau3 = 0.003;
tau4 = 0.09;
H_RI = (1+tau1*s)*exp(-delta*s)/((1+tau3*s)*(1+tau4*s));
Hd_RI = c2d(H_RI,1/Fs);
[num_RI_temp,den_RI_temp] = tfdata(Hd_RI);
num = cell2mat(num_RI_temp);
den = cell2mat(den_RI_temp);

num1 = num(1);
num2 = num(2);
num3 = num(3);
den1 = den(1);
den2 = den(2);
den3 = den(3);

FR_RI_temp = zeros(1,length(time));
FR_RI = zeros(1,length(time));
tic
for t = 1:length(time)
    U_mat(:,t) = U(t);
    %tic
    if t > 5
        %[FR_RI,FR_RI_temp] = RenshawOutput(FR_RI,FR_RI_temp,U,t,num,den);
        
        FR_RI_temp(t) = (num3.*U(t-2) + num2.*U(t-1) + num1.*U(t)...
            - den3.*FR_RI_temp(t-2) - den2.*FR_RI_temp(t-1))/den1;
        FR_RI(t) = FR_RI_temp(t);
    end
    %toc
    
end
toc
figure(1)
plot(time,U)
hold on 
plot(time,FR_RI_temp)
plot(time,FR_RI)

function [FR_RI,FR_RI_temp] = RenshawOutput(FR_RI,FR_RI_temp,ND,index,num,den)
num1 = num(1);
num2 = num(2);
num3 = num(3);
den1 = den(1);
den2 = den(2);
den3 = den(3);

FR_RI_temp(index) = (num3.*ND(index-2) + num2.*ND(index-1) + num1.*ND(index)...
    - den3.*FR_RI_temp(index-2) - den2.*FR_RI_temp(index-1))/den1;
FR_RI(index) = FR_RI_temp(index);
% if FR_RI(:,index) < 0
%     FR_RI(:,index) = 0;
% end

end