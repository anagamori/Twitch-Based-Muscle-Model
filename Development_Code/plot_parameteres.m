%close all
clc
clear all

code_folder = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Development_Code';
data_folder = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Model Parameters/Model_11';

N_MU = 200;
CT = zeros(N_MU,1);
t2t = zeros(N_MU,1);
FR_half = zeros(N_MU,1);
ratio = zeros(N_MU,1);
Af = zeros(N_MU,52);
fusion = zeros(N_MU,52);
parameter_Matrix = zeros(N_MU,14);
fusion_half = zeros(N_MU,1);
Af_f = zeros(N_MU,1);
fusion_f = zeros(N_MU,1);
f_exp = [2 4 6 8 10 12.5 14 16 18 20 22 25 28 30 33.3 40 50 66.6 80 100];
force_exp = [2.5723 7.717 16.399 29.582 44.373 50.482 62.379 63.344 67.203 73.234 77.492 83.923 87.46 88.103 91.318 96.141 97.428 100.32 100 98.392];
[~,loc] = min(abs(force_exp-50));
f_half_exp = f_exp(loc);
fusion_exp = [0 17.465 36.62 58.31 74.93 85.915 90.141 92.394 94.93 96.056 96.62 97.465 97.465 97.746 98.028 97.746 97.465 97.465 97.746 96.62];


for i = 1:N_MU
    
    cd('/Users/akira/Documents/Github/Twitch-Based-Muscle-Model/Development_Code/Data')
    load(['Data_v2_' num2str(i)])
    load(['phi_v2_' num2str(i)])
    cd(code_folder)
    
    CT(i) = Data{2,1};
    t2t(i) = Data{2,5};
    FR_half(i) = Data{2,6};
    ratio(i) = 1/FR_half(i)*1000/CT(i);
    
    Af(i,:) = Data{2,10}';
    fusion(i,:) = Data{2,11}';
    
    parameter = Data{2,12};
    parameter(13:14) = phi;
    parameter_Matrix(i,:) = parameter';
    
    figure(5)
    plot(Data{2,9},Data{2,10})
    xlim([0 3])
    hold on
    
    figure(3)
    plot(Data{2,9},Data{2,11}*100)
    xlim([0 3])
    hold on
    
    figure(4)
    ax_4 = plot(Data{2,10}*100,Data{2,11}*100,'color',[11,19,43]/255);
    ax_4.Color(4) = 0.5;
    hold on
    
    f_int = min(Data{2,9}):0.001:max(Data{2,9});
    fusion_int = spline(Data{2,9},Data{2,11}*100,f_int);
    [~,loc] = min(abs(f_int-1.1));
    fusion_f(i) = fusion_int(loc);
    
    Af_int = spline(Data{2,9},Data{2,10}*100,f_int);
    Af_f(i) = Af_int(loc);
    
    [~,loc_2] = min(abs(fusion_int-50));
    fusion_half(i) = f_int(loc_2);
end

figure(1)
plot(CT,t2t,'o')

figure(2)
plot(CT,1./FR_half*1000,'o')

figure(3)
plot(f_exp./f_half_exp,fusion_exp,'LineWidth',1,'color','k')

x_scale = 20/2.9;
y_scale =20/2.1;

x = [0.6 2.0 4.6 7 9.5 11.2 12.7 13.5 13.8];
y = [0.25 3 5 7.2 8.8 9.4 9.5 10.1 10.1];

Force = x*x_scale;
Fusion = y*y_scale;
%%
figure(4)
hold on
plot(Force,Fusion,'LineWidth',2,'color',[252,163,17]/255)
plot(force_exp,fusion_exp,'LineWidth',2,'color',[252,163,17]/255)
xlabel('Activation (%)')
ylabel('Fusion (%)')
set(gca,'TickDir','out');
set(gca,'box','off')

%%
figure(6)
plot(mean(Af)*100,mean(fusion)*100,'LineWidth',1,'color','k')
hold on
plot(Force,Fusion,'LineWidth',2,'color',[252,163,17]/255)

% %%
% cd('/Users/akira/Documents/Github/Twitch-Based-Muscle-Model/Development_Code/Data')
% save('CT','CT')
% save('t2t','t2t')
% save('FR_half','FR_half')
% %save('pool_parameter_matrix','parameter_Matrix')
% save('fusion','fusion')
% save('Af','Af')
% save('ratio','ratio')
% cd(code_folder)
