%==========================================================================
% excitation2dischargeRate.m
% Author: Akira Nagamori
% Last update: 4/21/19
% Descriptions:
%   Plot tiwtch amplitude of all motor units 
%   Used to generatee Fig XX in the manuscript
%==========================================================================
close all
cd('/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Model Parameters/Model_CTvsPTi');
%% Peak tension of muscle
density = 1.06; %
L0 = 6.8; % optimal muscle length [cm]
mass = 0.01; % muscle mass [kg]
PCSA = (mass*1000)/(density*L0); % PCSA of muscle
sigma = 31.8; % specific tension
F0 = PCSA * sigma; % maximal force

%% Number of motor unit
N_MU = 300; % number of motor units
i_MU = 1:N_MU; % index for motor units

%% Contraction time
% Generate a distribution of contraction time across motor units based on
% Rayleigh distribution
% rng(1)
% min_CT = 20; %minimum contraction time [ms]
% CT = round(raylrnd(23,1,N_MU)+min_CT); %contraction time of individual motor units [ms]
% CT_sorted = sort(CT,'descend');
load('CT_vec')
modelParameter.CT = CT_vec;

%% Peak tetanic force
RP_MU = 25; %range of peak tension across motor untis in unit of fold
b_MU = log(RP_MU)/N_MU; %coefficient to establish a range of twich force values
P_MU = exp(b_MU*i_MU); %force generated by a motor unit as a function of its recruitment threshold
PTi = P_MU./sum(P_MU)*F0; % peak tetanic force for individual units

%% Fractional PSCA
F_pcsa_slow = 0.3; % fractional PSCA of slow-twitch motor units (0-1)
[~, index_slow] = min(abs(cumsum(PTi) - F0*F_pcsa_slow)); rng(1)
%index_slow = 196;

%% Assign peak tetanic force into each unit
% shuffle indexes within each fiber type with respect to contraction time
% this will allow us to randomly assign peak tetanic tension to each motor
% unit with different contraction time
rng(1)
R_slow = randperm(index_slow);
index_fast = index_slow+1:N_MU;
R_fast_temp = randperm(length(index_fast));
R_fast = index_fast(R_fast_temp);
index_MU_PTi = [R_slow R_fast]; % vector of indexes to match peak tetanic tension to appropriate contraction time
PTi_new = PTi; %(index_MU_PTi);

%% Recruitment threshold
% Find recruitment threshold for individual units using exponential fit
% Recruitment threshold is correlated to peak tetanic tension
%   Use index_MU_PTi to appropriately index each MU
Ur = 0.5; % recruitment threshold for the lastly recruited motor unit
Ur_1 = 0.01; % reruitment threshold for the first unit
f_RT = fit([1 N_MU]',[Ur_1 Ur]','exp1');
coeffs_f_RT = coeffvalues(f_RT);
U_th = coeffs_f_RT(1)*exp(coeffs_f_RT(2)*i_MU); % the resulting recruitment threshold for individual units
U_th_new = U_th; %(index_MU_PTi);
[~,loc_max_U_th] = max(U_th_new);

%% FR_half for individual motor units
load('FR_half')
MDR = FR_half/2;
PDR = FR_half*2;

[~,index_DR_dif] = max(PDR-MDR);
%g_e = (PDR(index_DR_dif)-MDR(index_DR_dif))/(1-U_th_new(index_DR_dif));
g_e = (2-0.5)./(1-U_th_new(end));
%g_e = 115.1750;
%% 
load('t2t')
Pti = PTi_new.*t2t;
%%
cd('/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model')

%%
%randperm(300,50)
figure(1)
plot(CT_vec(1:196),PTi_new(1:196),'o','color',[36 123 160]/255,'LineWidth',1)
hold on 
plot(CT_vec(197:end),PTi_new(197:end),'o','color',[255 22 84]/255,'LineWidth',1)
xlabel('Contraction Time (ms)','FontSize',14)
ylabel('Peak Tetanic Force (N)','FontSize',14)
set(gca,'TickDir','out');
set(gca,'box','off')
ax = gca;
ax.FontSize = 10;

%%
mdl_CT_t2t_slow = fitlm(CT_vec(1:196),t2t(1:196));
betahat_CT_t2t_slow = mdl_CT_t2t_slow.Coefficients.Estimate;
mdl_CT_t2t_fast = fitlm(CT_vec(197:300),t2t(197:300));
betahat_CT_t2t_fast = mdl_CT_t2t_fast.Coefficients.Estimate;

figure(2)
plot(CT_vec(1:196),t2t(1:196),'o','color',[36 123 160]/255,'LineWidth',1)
hold on 
plot(CT_vec(1:196),betahat_CT_t2t_slow(1)+CT_vec(1:196)*betahat_CT_t2t_slow(2),'k','LineWidth',1)
plot(CT_vec(197:end),t2t(197:end),'o','color',[255 22 84]/255,'LineWidth',1)
plot(CT_vec(197:end),betahat_CT_t2t_fast(1)+CT_vec(197:end)*betahat_CT_t2t_fast(2),'k','LineWidth',1)
text(mean(CT_vec(1:196)),0.6,num2str(sqrt(mdl_CT_t2t_slow.Rsquared.Ordinary)))
text(mean(CT_vec(197:end)),0.6,num2str(sqrt(mdl_CT_t2t_fast.Rsquared.Ordinary)))
xlabel('Contraction Time (ms)','FontSize',14)
ylabel('Twitch-Tetanus Ratio','FontSize',14)
ylim([0 0.65])
set(gca,'TickDir','out');
set(gca,'box','off')
ax = gca;
ax.FontSize = 10;

%% 
figure(3)
histogram(t2t(1:196),0:0.05:0.6,'FaceColor',[36 123 160]/255)
hold on
histogram(t2t(197:end),0:0.05:0.6,'FaceColor',[255 22 84]/255)
xlabel('Twitch-Tetanus Ratio','FontSize',14)
ylabel('Counts','FontSize',14)
legend('Slow-twitch','Fast-twitch')
set(gca,'TickDir','out');
set(gca,'box','off')
ax = gca;
ax.FontSize = 10;

%%
figure(4)
plot(U_th_new(1:196),PTi_new(1:196),'o','color',[36 123 160]/255,'LineWidth',1)
hold on 
plot(U_th_new(197:end),PTi_new(197:end),'o','color',[255 22 84]/255,'LineWidth',1)
xlabel('Recruitment Threshold (Fraction of Maximum Activation)','FontSize',14)
ylabel('Peak Tetanic Force (N)','FontSize',14)
set(gca,'TickDir','out');
set(gca,'box','off')
ax = gca;
ax.FontSize = 10;

%%
mdl_CT_FR_half = fitlm(CT_vec,FR_half);
betahat_CT_FR_half = mdl_CT_FR_half.Coefficients.Estimate;

figure(5)
plot(CT_vec(1:196),FR_half(1:196),'o','color',[36 123 160]/255,'LineWidth',1)
hold on 
plot(CT_vec(197:end),FR_half(197:end),'o','color',[255 22 84]/255,'LineWidth',1)
plot(CT_vec,betahat_CT_FR_half(1)+CT_vec*betahat_CT_FR_half(2),'k','LineWidth',1)
text(mean(CT_vec),45,num2str(sqrt(mdl_CT_FR_half.Rsquared.Ordinary)))
xlabel('Contraction Time (ms)','FontSize',14)
ylabel('f_{0.5} (Hz)','FontSize',14)
set(gca,'TickDir','out');
set(gca,'box','off')
ax = gca;
ax.FontSize = 10;

%%
figure(6)
histogram(PTi(1:196),0:0.02:0.5,'FaceColor',[36 123 160]/255)
hold on
histogram(PTi(197:end),0:0.02:0.5,'FaceColor',[255 22 84]/255)
xlabel('Peak Tetanic Force','FontSize',14)
ylabel('Counts','FontSize',14)
legend('Slow-twitch','Fast-twitch')
set(gca,'TickDir','out');
set(gca,'box','off')
ax = gca;
ax.FontSize = 10;

%%
figure(7)
histogram(CT_vec(1:196),10:10:120,'FaceColor',[36 123 160]/255)
hold on
histogram(CT_vec(197:end),10:10:120,'FaceColor',[255 22 84]/255)
xlabel('Contraction Time (s)','FontSize',14)
ylabel('Counts','FontSize',14)
legend('Slow-twitch','Fast-twitch')
set(gca,'TickDir','out');
set(gca,'box','off')
ax = gca;
ax.FontSize = 10;

