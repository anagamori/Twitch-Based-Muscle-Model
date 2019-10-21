%==========================================================================
% dischargeRate2excitation_nonLinear.m
% Author: Akira Nagamori
% Last update: 6/27/19
% Descriptions:
%   Plot the relationship between excitation and discharge rate of
%   individual motor units
%   Used for c of Summary_Recruitment.pdf
%==========================================================================

cd('/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Model Parameters/Model_8');
load('modelParameter')
cd('/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Code for Figures')

%%
N_MU = modelParameter.N_MU;
MDR = modelParameter.MDR;
PDR = modelParameter.PDR;
U_th = modelParameter.U_th;
%%
g_e = modelParameter.g_e;
lamda =  modelParameter.lamda;
k_e =  modelParameter.k_e;
index_saturation = modelParameter.index_saturation;
U_th_t =  modelParameter.U_th_t;

%% Discharge rate of motor unit
U_vec = 0:0.001:1;
DR_mat = zeros(N_MU,length(U_vec));
DR_temp = zeros(N_MU,1);
for i = 1:length(U_vec)
    DR_MU = g_e.*(U_vec(i)-U_th)+MDR;   
%    DR_temp = MDR + lamda.*k_e.*(U_vec(i)-U_th_new);
    for n = 1:length(index_saturation)
        index = index_saturation(n);
        if U_vec(i) <= U_th_t(index)
            DR_temp(index) = MDR(index) + lamda(index).*k_e(index).*(U_vec(i)-U_th(index));
        else
            DR_temp(index) = PDR(index)-k_e(index)*(1-U_vec(i));
        end
    end
    
    DR_MU(index_saturation) = DR_temp(index_saturation);
    % DR_MU(1:index_slow) = DR_temp(1:index_slow);
    DR_MU(DR_MU<MDR) = NaN;
    DR_MU(DR_MU>PDR) = PDR(DR_MU>PDR);
    DR_mat(:,i) = DR_MU;
end

%%
%close all
%index_plot = 1:1:300; %[1 50 100 150 200 250 300];
index_plot = [2 10 38 40 60 80 83 100 111 120];
figure(1)
plot(U_vec*100,DR_mat(index_plot,:),'k','LineWidth',1)
xlabel('Activation (%Maximum)','FontSize',8)
ylabel('Discharge Rate (Hz)','FontSize',8)
set(gca,'TickDir','out');
set(gca,'box','off')
ax = gca;
ax.FontSize = 6;
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 3.34 3.34];
hold on
