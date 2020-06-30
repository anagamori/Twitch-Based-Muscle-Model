%==========================================================================
%PDRvsUth_nonlinear.m
% Author: Akira Nagamori
% Last update: 6/27/20
% Descriptions:
%   Plot the relationship between excitation and discharge rate of
%   individual motor units
%   Used to generate figure (onion_skinVsAHP)
%==========================================================================
close all
clear all
clc

cd('/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Development_Code/Data');
load('modelParameter_v3')
cd('/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Code for Figures')

%%
MDR = modelParameter.MDR;
PDR = modelParameter.PDR;
U_th = modelParameter.U_th;
%%
g_e = modelParameter.g_e;
lamda =  modelParameter.lamda;
k_e =  modelParameter.k_e;
index_saturation = modelParameter.index_saturation;
U_th_t =  modelParameter.U_th_t;
%%
comb = nchoosek(1:modelParameter.N_MU,2);
rng shuffle
for i = 1:1000
    index = randperm(size(comb,1),1);
    if U_th(comb(index,2)) >= U_th(comb(index,1))
        U_th_diff(i) = U_th(comb(index,2)) - U_th(comb(index,1));
        PDR_diff(i) = PDR(comb(index,2)) - PDR(comb(index,1));
        MDR_diff(i) = MDR(comb(index,2)) - MDR(comb(index,1));
    else
        U_th_diff(i) = U_th(comb(index,1)) - U_th(comb(index,2));
        PDR_diff(i) = PDR(comb(index,1)) - PDR(comb(index,2));
        MDR_diff(i) = MDR(comb(index,1)) - MDR(comb(index,2));
    end
    
end

%%
% length(find(PDR_diff>0))/length(PDR_diff)
% figure(1)
% scatterhist(U_th_diff*100,PDR_diff,'Color','k','MarkerSize',10,'Marker','.')
% xlim([-5 55])
% hold on
% plot([-5 55],[0 0],'r','linewidth',2)
% xlabel('Difference in Recruitment Threshold (%Maximum)','FontSize',8)
% ylabel('Difference in Peak Discharge Rate (Hz)','FontSize',8)
% set(gca,'TickDir','out');
% set(gca,'box','off')
% ax = gca;
% ax.FontSize = 6;
% fig = gcf;
% fig.PaperUnits = 'inches';
% fig.PaperPosition = [0 0 3.34 3.34];

%%
U_th_vec = U_th;
U_th_vec =  U_th_vec*100;
X_1 = [ones(length(U_th_vec),1) U_th_vec];
b_1 = X_1\PDR;
PDR_Calc = X_1*b_1;
Rsq_1 = 1 - sum((PDR - PDR_Calc).^2)/sum((PDR - mean(PDR)).^2)
[R_1,P_1] = corrcoef(U_th_vec,PDR)

figure(2)
scatter(U_th*100,PDR,'filled','k','LineWidth',0.5)
hold on
plot(U_th*100,PDR_Calc,'r','LineWidth',1)
xlim([-5 85])
xticks([0 0.1 0.2 0.3 0.4 0.5]*100)
xlabel('Recruitment Threshold (%Maximum)','FontSize',8)
ylabel('Discharge Rate (Hz)','FontSize',8)
set(gca,'TickDir','out');
set(gca,'box','off')
ax = gca;
ax.FontSize = 6;
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 3.34 3.34];
%%
clear U_th_diff
DR_temp = zeros(modelParameter.N_MU,1);
U = 0.2;
DR_MU = g_e.*(U-U_th)+MDR;

%    DR_temp = MDR + lamda.*k_e.*(U_vec(i)-U_th_new);
for n = 1:length(index_saturation)
    index = index_saturation(n);
    if U <= U_th_t(index)
        DR_temp(index) = MDR(index) + lamda(index).*k_e(index).*(U-U_th(index));
    else
        DR_temp(index) = PDR(index)-k_e(index)*(1-U);
    end
end

DR_MU(index_saturation) = DR_temp(index_saturation);
% DR_MU(1:index_slow) = DR_temp(1:index_slow);
DR_MU(DR_MU<MDR) = 0;
DR_MU(DR_MU>PDR) = PDR(DR_MU>PDR);

DR_vec = DR_MU(DR_MU>0);
U_th_vec = U_th(DR_MU>0);
U_th_vec = U_th_vec*100;

X = [ones(length(U_th_vec),1) U_th_vec];
b1 = X\DR_vec;
DR_Calc = X*b1;
Rsq = 1 - sum((DR_vec - DR_Calc).^2)/sum((DR_vec - mean(DR_vec)).^2)
[R,P] = corrcoef(U_th(DR_MU>0),DR_MU(DR_MU>0))

figure(3)
scatter(U_th_vec,DR_vec,'filled','k','LineWidth',0.5)
hold on
plot(U_th_vec,DR_Calc,'r','LineWidth',1)
xlim([-1 11])
xticks([0 0.05 0.1 0.15 0.2]*100)
xlabel('Recruitment Threshold (%Maximum)','FontSize',8)
ylabel('Discharge Rate (Hz)','FontSize',8)
set(gca,'TickDir','out');
set(gca,'box','off')
ax = gca;
ax.FontSize = 6;
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 3.34 3.34];

%%
% comb = nchoosek(find(DR_MU>0),2);
% rng shuffle
% for i = 1:1000
%     index = randperm(size(comb,1),1);
%     if U_th(comb(index,2)) >= U_th(comb(index,1))
%         U_th_diff(i) = U_th(comb(index,2)) - U_th(comb(index,1));
%         DR_diff(i) = DR_MU(comb(index,2)) - DR_MU(comb(index,1));
%         
%     else
%         U_th_diff(i) = U_th(comb(index,1)) - U_th(comb(index,2));
%         DR_diff(i) = DR_MU(comb(index,1)) - DR_MU(comb(index,2));
%     end
% end
% 
% figure(4)
% scatterhist(U_th_diff,DR_diff,'Color','k')
% hold on
% plot([-0.05 0.5],[0 0],'r','linewidth',2)
% xlabel('Difference in Recruitment Threshold (%Maximum)','FontSize',8)
% ylabel('Difference in Discharge Rate (Hz)','FontSize',8)
% set(gca,'TickDir','out');
% set(gca,'box','off')
% ax = gca;
% ax.FontSize = 6;
% fig = gcf;
% fig.PaperUnits = 'inches';
% fig.PaperPosition = [0 0 3.34 3.34];
% 
% length(find(DR_diff>0))/length(DR_diff)




