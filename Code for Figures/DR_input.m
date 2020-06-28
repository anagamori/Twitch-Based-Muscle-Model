%==========================================================================
%PDRvsUth_nonlinear.m
% Author: Akira Nagamori
% Last update: 6/27/19
% Descriptions:
%   Plot the relationship between excitation and discharge rate of
%   individual motor units
%   Used to generate figure (onion_skinVsAHP)
%==========================================================================
close all
clear all
clc

cd('/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Development_Code/Data');
load('modelParameter_v2')
cd('/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Code for Figures')

%%
x = [11.025  21.277 50.097 72.34 93.617];
DR_Connelly = [12.269 15.231 22.212 28.135 41.779];
error_Connelly = [14.596 19.356 28.452 35.856 50.029];
error_Connelly = error_Connelly - DR_Connelly;
	
%%
MDR = modelParameter.MDR;
PDR = modelParameter.PDR;
U_th = modelParameter.U_th;
index_slow = modelParameter.index_slow;

%%
g_e = modelParameter.g_e;
lamda =  modelParameter.lamda;
k_e =  modelParameter.k_e;
index_saturation = modelParameter.index_saturation;
U_th_t =  modelParameter.U_th_t;

U_vec = x/100; %[0.1 0.25 0.5 0.75 1];
DR_MU_all = zeros(200,length(U_vec));

for i = 1:length(U_vec)

DR_temp = zeros(modelParameter.N_MU,1);
U = U_vec(i);
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

DR_MU(DR_MU<MDR) = nan;
DR_MU(DR_MU>PDR) = PDR(DR_MU>PDR);

DR_MU_all(:,i) = DR_MU;

end



%%
figure(1)
% boxplot(DR_MU_all,'positions',[10 25 50 75 100])
% hold on 
% scatter(10*ones(size(DR_MU_all(:,1))).*(1+(rand(size(DR_MU_all(:,1)))-0.5)/10),DR_MU_all(:,1),'r','filled')

errorbar(x+1,nanmean(DR_MU_all),nanstd(DR_MU_all),'LineWidth',1,'Color',[37  65 178]/255)
hold on
errorbar(x,DR_Connelly,error_Connelly,'LineWidth',1,'Color','k')
xlabel('Synaptic Input (%Maximum)','FontSize',12)
ylabel('Dsicharge Rate (Hz)','FontSize',12)
legend('New Model','Connelly et al. 1999')
set(gca,'TickDir','out');
set(gca,'box','off')
% scatter(U_th_vec,DR_vec,'filled','k','LineWidth',0.5)
% hold on
% plot(U_th_vec,DR_Calc,'r','LineWidth',1)
% xlim([-1 11])
% xticks([0 0.05 0.1 0.15 0.2]*100)
% xlabel('Recruitment Threshold (%Maximum)','FontSize',8)
% ylabel('Discharge Rate (Hz)','FontSize',8)
% set(gca,'TickDir','out');
% set(gca,'box','off')
% ax = gca;
% ax.FontSize = 6;
% fig = gcf;
% fig.PaperUnits = 'inches';
% fig.PaperPosition = [0 0 3.34 3.34];




