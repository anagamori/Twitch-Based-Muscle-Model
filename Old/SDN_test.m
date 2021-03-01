%==========================================================================
% plotData_withTendon.m
% Author: Akira Nagamori
% Last update: 3/11/19
% Descriptions:
%==========================================================================
close all
clear all
clc

%%
condition = 'Model_6_20_CoV_50_Ur_Rec_3';
data_folder = ['/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Data/withTendon/' condition];
code_folder = '/Users/akiranagamori/Documents/Github/Twitch-Based-Muscle-Model';
figure_folder = '/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Figures';

Fs = 10000;
amp_vec = [0.05 0.1:0.1:1];
time =0:1/Fs:15;
mean_Force = zeros(10,length(amp_vec));
std_Force = zeros(10,length(amp_vec));
cov_Force = zeros(10,length(amp_vec));
var_Force = zeros(10,length(amp_vec));



trial_vec = 0:10; %0:10; %[0:6 8:10];
for k = 1:length(trial_vec)
    j = trial_vec(k);
    cd(data_folder)
    load(['Force_mat_' num2str(j) ' 2'])
    cd(code_folder)
    for i = 1:10
        Force =  Force_mat(i,:);
        mean_Force(i,j+1) = mean(Force(5*Fs+1:end));
        std_Force(i,j+1) = std(Force(5*Fs+1:end));
        cov_Force(i,j+1) =  std_Force(i,j+1)/mean_Force(i,j+1)*100;
        var_Force(i,j+1) =  var(Force(5*Fs+1:end));
        
    end
    
    figure(11)
    plot(time,Force_mat)
    hold on
        
    %clear Force_mat
end

%%
amp_vec = amp_vec.^2*100;
mean_mean_Force = mean(mean_Force);
MVC = mean_mean_Force(end);
var_Force = var_Force./MVC*100;
%%
p = polyfit(amp_vec(1:6),var_Force(1:6),1);


figure(2)
errorbar(amp_vec,mean(var_Force),std(var_Force),'LineWidth',2,'Color','b');
hold on 
plot(amp_vec,p(1)+p(2)*amp_vec,'--k')
xlabel('Activation Level/Mean Force (%)','FontSize',14)
ylabel('Variance (%MVC)','FontSize',14)

%%
c = sqrt(p(2))