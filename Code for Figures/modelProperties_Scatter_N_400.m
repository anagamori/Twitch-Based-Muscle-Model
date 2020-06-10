%==========================================================================
% modelProperties_Scatter_v2.m
% Author: Akira Nagamori
% Last update: 5/26/2020
% Descriptions:
%   Plot tiwtch amplitude of all motor units 
%   Used to generate Model_Characteristics in the manuscript
%==========================================================================
close all
clear all
clc

cd('/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Model Parameters/Model_11_N_400');
load('modelParameter')
load('CT_vec_new')
load('FR_half')
load('t2t')
load('index_MU')
cd('/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Code for Figures')

%%
N_MU = modelParameter.N_MU;

CT = CT_vec';
%t2t = t2t(index_MU)';
FR_half = modelParameter.FR_half;
MDR = FR_half/2;
PDR = FR_half*2;
PTi = modelParameter.PTi;
index_slow = modelParameter.index_slow;
U_th = modelParameter.U_th;

%%
%randperm(300,50)

figure(1)
plot(CT(1:index_slow),PTi(1:index_slow),'o','color',[36 123 160]/255,'LineWidth',1)
hold on 
plot(CT(index_slow+1:end),PTi(index_slow+1:end),'o','color',[255 22 84]/255,'LineWidth',1)
xlabel('Contraction Time (ms)','FontSize',14)
ylabel('Peak Tetanic Force (N)','FontSize',14)
set(gca,'TickDir','out');
set(gca,'box','off')
ax = gca;
ax.FontSize = 10;

%%
% mdl_CT_t2t_slow = fitlm(CT(1:index_slow),t2t(1:index_slow));
% betahat_CT_t2t_slow = mdl_CT_t2t_slow.Coefficients.Estimate;
% mdl_CT_t2t_fast = fitlm(CT(index_slow+1:N_MU),t2t(index_slow+1:N_MU));
% betahat_CT_t2t_fast = mdl_CT_t2t_fast.Coefficients.Estimate;
% 
% figure(2)
% plot(CT(1:index_slow),t2t(1:index_slow),'o','color',[36 123 160]/255,'LineWidth',1)
% hold on 
% plot(CT(1:index_slow),betahat_CT_t2t_slow(1)+CT(1:index_slow)*betahat_CT_t2t_slow(2),'k','LineWidth',1)
% plot(CT(index_slow+1:end),t2t(index_slow+1:end),'o','color',[255 22 84]/255,'LineWidth',1)
% plot(CT(index_slow+1:end),betahat_CT_t2t_fast(1)+CT(index_slow+1:end)*betahat_CT_t2t_fast(2),'k','LineWidth',1)
% text(mean(CT(1:index_slow)),0.6,num2str(sqrt(mdl_CT_t2t_slow.Rsquared.Ordinary)))
% text(mean(CT(index_slow+1:end)),0.6,num2str(sqrt(mdl_CT_t2t_fast.Rsquared.Ordinary)))
% xlabel('Contraction Time (ms)','FontSize',14)
% ylabel('Twitch-Tetanus Ratio','FontSize',14)
% ylim([0 0.65])
% set(gca,'TickDir','out');
% set(gca,'box','off')
% ax = gca;
% ax.FontSize = 10;

%% 
% figure(3)
% histogram(t2t(1:index_slow),0:0.05:0.6,'FaceColor',[36 123 160]/255)
% hold on
% histogram(t2t(index_slow+1:end),0:0.05:0.6,'FaceColor',[255 22 84]/255)
% xlabel('Twitch-Tetanus Ratio','FontSize',14)
% ylabel('Counts','FontSize',14)
% legend('Slow-twitch','Fast-twitch')
% set(gca,'TickDir','out');
% set(gca,'box','off')
% ax = gca;
% ax.FontSize = 10;
% 
%%
figure(4)
plot(U_th(1:index_slow),PTi(1:index_slow),'o','color',[36 123 160]/255,'LineWidth',1)
hold on 
plot(U_th(index_slow+1:end),PTi(index_slow+1:end),'o','color',[255 22 84]/255,'LineWidth',1)
xlabel('Recruitment Threshold (Fraction of Maximum Activation)','FontSize',14)
ylabel('Peak Tetanic Force (N)','FontSize',14)
set(gca,'TickDir','out');
set(gca,'box','off')
ax = gca;
ax.FontSize = 10;

%%
ISI = 1./FR_half*1000;
mdl_CT_FR_half = fitlm(CT,ISI);
betahat_CT_FR_half = mdl_CT_FR_half.Coefficients.Estimate;

figure(5)
plot(CT(1:index_slow),ISI(1:index_slow),'o','color',[36 123 160]/255,'LineWidth',1)
hold on 
plot(CT(index_slow+1:end),ISI(index_slow+1:end),'o','color',[255 22 84]/255,'LineWidth',1)
plot(CT,betahat_CT_FR_half(1)+CT*betahat_CT_FR_half(2),'k','LineWidth',1)
text(80,90,num2str(sqrt(mdl_CT_FR_half.Rsquared.Ordinary)))
text(80,80,num2str(betahat_CT_FR_half(2)))
%ylim([0 45])
xlabel('Contraction Time (ms)','FontSize',14)
ylabel('Stimulus Interval at f_{0.5} (ms)','FontSize',14)
set(gca,'TickDir','out');
set(gca,'box','off')
ax = gca;
ax.FontSize = 10;

%%
figure(6)
histogram(PTi(1:index_slow),0:1:12,'FaceColor',[36 123 160]/255)
hold on
histogram(PTi(index_slow+1:end),0:1:12,'FaceColor',[255 22 84]/255)
xlabel('Peak Tetanic Force','FontSize',14)
ylabel('Counts','FontSize',14)
legend('Slow-twitch','Fast-twitch')
set(gca,'TickDir','out');
set(gca,'box','off')
ax = gca;
ax.FontSize = 10;

%%
x = 22.5:5:87.5;
frac_Cuseum = [0.016941 0.056471 0.13882 0.19859 0.17129 0.16094 0.11059 0.049882 0.024471 0.020706 0.024471 0.016941 0.0056471 0.0089412];

int = 20:5:90;
count_slow = zeros(1,length(x));
count_fast = zeros(1,length(x));
for i = 1:length(int)-1
    count_slow(i) = length(find(CT(1:index_slow)>=int(i)&CT(1:index_slow)<int(i+1)));
    count_fast(i) = length(find(CT(index_slow+1:end)>=int(i)&CT(index_slow+1:end)<int(i+1)));
end
frac_slow = count_slow./(sum(count_slow)+sum(count_fast));
frac_fast = count_fast./(sum(count_slow)+sum(count_fast));
figure(7)
h = bar(x,frac_slow,'hist');
set(h,'FaceColor',[36 123 160]/255,'FaceAlpha',0.5);
hold on 
h2 = bar(x,frac_fast,'hist');
set(h2,'FaceColor',[255 22 84]/255,'FaceAlpha',0.5);
plot(x,frac_Cuseum,'k','LineWidth',1)
plot(x,frac_Cuseum,'o', 'Color','k','MarkerFaceColor', 'k')
xlabel('Contraction Time (s)','FontSize',10)
ylabel('Proportion','FontSize',14)
legend('Slow-twitch','Fast-twitch','Cutsem et al. 1993')
set(gca,'TickDir','out');
set(gca,'box','off')
ax = gca;
ax.FontSize = 10;

%%
x = 5:2:15;
frac_Cutsem = [0.17846 0.20369 0.23015 0.17908 0.11508 0.038769];
int = 2:2:20;
count_slow = zeros(1,length(int)-1);
count_fast = zeros(1,length(int)-1);
for i = 1:length(int)-1
    count_slow(i) = length(find(MDR(1:index_slow)>=int(i)&MDR(1:index_slow)<int(i+1)));
    count_fast(i) = length(find(MDR(index_slow+1:end)>=int(i)&MDR(index_slow+1:end)<int(i+1)));
end
frac_slow = count_slow./(sum(count_slow)+sum(count_fast));
frac_fast = count_fast./(sum(count_slow)+sum(count_fast));

figure(8)
h = bar(3:2:19,frac_slow,'hist');
set(h,'FaceColor',[36 123 160]/255,'FaceAlpha',0.5);
set(gca,'XTick',2:2:20)
hold on 
h2 = bar(3:2:19,frac_fast,'hist');
set(h2,'FaceColor',[255 22 84]/255,'FaceAlpha',0.5);
plot(x,frac_Cutsem,'k','LineWidth',1)
plot(x,frac_Cutsem,'o', 'Color','k','MarkerFaceColor', 'k')
xlabel('MDR (Hz)','FontSize',10)
ylabel('Proportion','FontSize',14)
legend('Slow-twitch','Fast-twitch','Cutsem et al. 1993')
set(gca,'TickDir','out');
set(gca,'box','off')
ax = gca;
ax.FontSize = 10;

%%
x = 15:10:65;
frac_Cutsem = [0.071464 0.26427 0.30149 0.14516 0.12134 0.026055];

int = 10:10:80;
count_slow = zeros(1,length(int)-1);
count_fast = zeros(1,length(int)-1);
for i = 1:length(int)-1
    count_slow(i) = length(find(PDR(1:index_slow)>=int(i)&PDR(1:index_slow)<int(i+1)));
    count_fast(i) = length(find(PDR(index_slow+1:end)>=int(i)&PDR(index_slow+1:end)<int(i+1)));
end
frac_slow = count_slow./(sum(count_slow)+sum(count_fast));
frac_fast = count_fast./(sum(count_slow)+sum(count_fast));

figure(9)
h = bar(15:10:75,frac_slow,'hist');
set(h,'FaceColor',[36 123 160]/255,'FaceAlpha',0.5);
set(gca,'XTick',10:10:80)
hold on 
h2 = bar(15:10:75,frac_fast,'hist');
set(h2,'FaceColor',[255 22 84]/255,'FaceAlpha',0.5);
plot(x,frac_Cutsem,'k','LineWidth',1)
plot(x,frac_Cutsem,'o', 'Color','k','MarkerFaceColor', 'k')
xlabel('PDR (Hz)','FontSize',10)
ylabel('Proportion','FontSize',14)
legend('Slow-twitch','Fast-twitch','Cutsem et al. 1993')
set(gca,'TickDir','out');
set(gca,'box','off')
ax = gca;
ax.FontSize = 10;

