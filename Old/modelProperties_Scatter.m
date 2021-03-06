%==========================================================================
% excitation2dischargeRate.m
% Author: Akira Nagamori
% Last update: 4/21/19
% Descriptions:
%   Plot tiwtch amplitude of all motor units 
%   Used to generate Model_Characteristics in the manuscript
%==========================================================================
close all
clear all
clc

cd('/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Model Parameters/Model_8');
load('modelParameter')
load('CT_vec')
load('t2t')
cd('/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Code for Figures')

%%
N_MU = modelParameter.N_MU;
index_MU = [1:300/100:300 300:-3*5:2];
index_MU = sort(index_MU);


CT = CT_vec(index_MU)';
t2t = t2t(index_MU)';
FR_half = modelParameter.FR_half;
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
mdl_CT_t2t_slow = fitlm(CT(1:index_slow),t2t(1:index_slow));
betahat_CT_t2t_slow = mdl_CT_t2t_slow.Coefficients.Estimate;
mdl_CT_t2t_fast = fitlm(CT(index_slow+1:N_MU),t2t(index_slow+1:N_MU));
betahat_CT_t2t_fast = mdl_CT_t2t_fast.Coefficients.Estimate;

figure(2)
plot(CT(1:index_slow),t2t(1:index_slow),'o','color',[36 123 160]/255,'LineWidth',1)
hold on 
plot(CT(1:index_slow),betahat_CT_t2t_slow(1)+CT(1:index_slow)*betahat_CT_t2t_slow(2),'k','LineWidth',1)
plot(CT(index_slow+1:end),t2t(index_slow+1:end),'o','color',[255 22 84]/255,'LineWidth',1)
plot(CT(index_slow+1:end),betahat_CT_t2t_fast(1)+CT(index_slow+1:end)*betahat_CT_t2t_fast(2),'k','LineWidth',1)
text(mean(CT(1:index_slow)),0.6,num2str(sqrt(mdl_CT_t2t_slow.Rsquared.Ordinary)))
text(mean(CT(index_slow+1:end)),0.6,num2str(sqrt(mdl_CT_t2t_fast.Rsquared.Ordinary)))
xlabel('Contraction Time (ms)','FontSize',14)
ylabel('Twitch-Tetanus Ratio','FontSize',14)
ylim([0 0.65])
set(gca,'TickDir','out');
set(gca,'box','off')
ax = gca;
ax.FontSize = 10;

%% 
figure(3)
histogram(t2t(1:index_slow),0:0.05:0.6,'FaceColor',[36 123 160]/255)
hold on
histogram(t2t(index_slow+1:end),0:0.05:0.6,'FaceColor',[255 22 84]/255)
xlabel('Twitch-Tetanus Ratio','FontSize',14)
ylabel('Counts','FontSize',14)
legend('Slow-twitch','Fast-twitch')
set(gca,'TickDir','out');
set(gca,'box','off')
ax = gca;
ax.FontSize = 10;

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
mdl_CT_FR_half = fitlm(CT,FR_half);
betahat_CT_FR_half = mdl_CT_FR_half.Coefficients.Estimate;

figure(5)
plot(CT(1:index_slow),FR_half(1:index_slow),'o','color',[36 123 160]/255,'LineWidth',1)
hold on 
plot(CT(index_slow+1:end),FR_half(index_slow+1:end),'o','color',[255 22 84]/255,'LineWidth',1)
plot(CT,betahat_CT_FR_half(1)+CT*betahat_CT_FR_half(2),'k','LineWidth',1)
text(mean(CT),45,num2str(sqrt(mdl_CT_FR_half.Rsquared.Ordinary)))
ylim([0 45])
xlabel('Contraction Time (ms)','FontSize',14)
ylabel('f_{0.5} (Hz)','FontSize',14)
set(gca,'TickDir','out');
set(gca,'box','off')
ax = gca;
ax.FontSize = 10;

%%
figure(6)
histogram(PTi(1:index_slow),0:0.1:2,'FaceColor',[36 123 160]/255)
hold on
histogram(PTi(index_slow+1:end),0:0.1:2,'FaceColor',[255 22 84]/255)
xlabel('Peak Tetanic Force','FontSize',14)
ylabel('Counts','FontSize',14)
legend('Slow-twitch','Fast-twitch')
set(gca,'TickDir','out');
set(gca,'box','off')
ax = gca;
ax.FontSize = 10;

%%
figure(7)
histogram(CT(1:index_slow),10:10:120,'FaceColor',[36 123 160]/255)
hold on
histogram(CT(index_slow+1:end),10:10:120,'FaceColor',[255 22 84]/255)
xlabel('Contraction Time (s)','FontSize',14)
ylabel('Counts','FontSize',14)
legend('Slow-twitch','Fast-twitch')
set(gca,'TickDir','out');
set(gca,'box','off')
ax = gca;
ax.FontSize = 10;

