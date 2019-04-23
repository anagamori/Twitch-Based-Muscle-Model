close all
clear all
clc

%%
code_folder = '/Users/akiranagamori/Documents/Github/Twitch-Based-Muscle-Model';
figure_folder = '/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Figures';

condition = '10_CoV_50_Ur_Rec_2_CTvsPTi';
Fs = 2000;
time =0:1/Fs:15;
data_folder = ['/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Data/withTendon/' condition];
cd(data_folder)
load('mean_Force')
load('std_Force')
load('cov_Force')
load('mean_pxx')
cd(code_folder)
color_code = [215 25 28]/255;

cd('/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Model Parameters/Model_CTvsPTi')
load('modelParameter')
cd(code_folder)

Force = [0 mean(mean_Force)];
Force = Force./Force(end);
act = 0:0.1:1;
act_int = 0:0.0001:1;
Force_int = interp1(act,Force,act_int);

for i = 1:length(modelParameter.U_th_new)
    [loc,index] = min(abs(modelParameter.U_th_new(i) - act_int));
    U_th_force(i) = Force_int(index);
end

figure(1)
plot(act_int,Force_int)

figure(2)
plot(modelParameter.U_th_new,'o')
hold on
plot(U_th_force,'o')