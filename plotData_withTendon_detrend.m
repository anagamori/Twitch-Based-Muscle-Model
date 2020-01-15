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
condition = 'Model_8_20_CoV_50_Ur_Rec_3';
data_folder = ['/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Data/withTendon/' condition];
code_folder = '/Users/akira/Documents/Github/Twitch-Based-Muscle-Model';
figure_folder = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Figures';

 Fs = 10000;
        time = 0:1/Fs:15;
amp_vec = [0.025 0.05 0.1:0.1:1];
nTrial = 10;
mean_Force = zeros(nTrial,length(amp_vec));
std_Force = zeros(nTrial,length(amp_vec));
cov_Force = zeros(nTrial,length(amp_vec));
cov_Force_dt = zeros(nTrial,length(amp_vec));
pxx = zeros(nTrial,201);
mean_pxx = zeros(length(amp_vec),201);

trial_vec = 0:10; %-1:10; %0:10; %[0:6 8:10];
for k = 1:length(trial_vec)
    j = trial_vec(k);
   
    cd(data_folder)
    load(['Force_mat_' num2str(j)])
    cd(code_folder)
    for i = 1:10
        Force_temp =  Force_mat(i,:);
        Force = Force_temp(5*Fs+1:end);
        Force_dt = detrend(Force,1,1*Fs);
        bp_time = 1:1*Fs:length(Force_dt);
        cov_Force_segment = zeros(1,length(bp_time)-1);
        for n = 1:length(bp_time)-1
            Force_segment = Force(bp_time(n):bp_time(n+1));
            Force_segment_dt = detrend(Force_segment,1);
            mean_Force_segment = mean(Force_segment);
            std_Force_segment = std(Force_segment_dt);
            cov_Force_segment(n) = std_Force_segment/mean_Force_segment*100;
        end
        mean_Force(i,j+2) = mean(Force);
        std_Force(i,j+2) = std(Force);
        cov_Force(i,j+2) =  std_Force(i,j+2)/mean_Force(i,j+2)*100;
        cov_Force_dt(i,j+2) =  mean(cov_Force_segment);
              
        [pxx(i,:),f] = pwelch(Force-mean(Force),rectwin(length(Force)),0,0:0.5:100,Fs,'power');
    end
    
    figure(11)
    plot(time,Force_mat)
    hold on
        
    mean_pxx(j+2,:) = mean(pxx);
    %clear Force_mat
end

cd(data_folder)
save('mean_Force','mean_Force')
save('std_Force','std_Force')
save('cov_Force','cov_Force')
save('cov_Force_dt','cov_Force_dt')
save('mean_pxx','mean_pxx')
cd(code_folder)
%%
%close all
mean_mean_Force = mean(mean_Force);
figure(1)
plot([0 amp_vec],[0 mean(mean_Force)]./mean_mean_Force(end),'LineWidth',2)
xlabel('Activation (%)','FontSize',14)
ylabel('Force (%MVC)','FontSize',14)
set(gca,'TickDir','out');
set(gca,'box','off')
% cd (figure_folder)
% saveas(gcf,['activation2meanForce_withTendon_' condition],'pdf')
% cd (code_folder)

%%
figure(2)
errorbar(amp_vec,mean(std_Force),std(std_Force),'LineWidth',2,'Color','b');
xlabel('Activation Level/Mean Force (%)','FontSize',14)
ylabel('SD (N)','FontSize',14)
hold on 
errorbar(mean(mean_Force)./mean_mean_Force(end),mean(std_Force),std(std_Force),'LineWidth',2,'Color','r');
legend('Activation Level','Force Level','Location','northwest')
yticks([0.05 0.1 0.15 0.2 0.25])
set(gca,'TickDir','out');
set(gca,'box','off')
% cd (figure_folder)
% saveas(gcf,['meanForce2SD_wtihTendon_' condition],'pdf')
% cd (code_folder)

%%
figure(3)
errorbar(amp_vec,mean(cov_Force),std(cov_Force),'LineWidth',2,'Color','b');
xlabel('Activation Level/Mean Force (%)','FontSize',14)
ylabel('CoV (%)','FontSize',14)
hold on 
errorbar(mean(mean_Force)./mean_mean_Force(end),mean(cov_Force),std(cov_Force),'LineWidth',2,'Color','r');
legend('Activation Level','Force Level')
set(gca,'TickDir','out');
set(gca,'box','off')
% cd (figure_folder)
% saveas(gcf,['meanForce2CoV_withTendon_' condition],'pdf')
% cd (code_folder)

%%
figure(4)
plot(f,mean_pxx([1 3 5 8],:),'LineWidth',2)
xlabel('Frequency (Hz)','FontSize',14)
ylabel('Power (N^2)','FontSize',14)
set(gca,'TickDir','out');
set(gca,'box','off')
legend('10% MVC','30% MVC','50% MVC','80% MVC')
% cd (figure_folder)
% saveas(gcf,['pxx_withTendon_' condition],'pdf')
% cd (code_folder)