close all
clc
clear all

code_folder = '/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Development_Code';
data_folder = '/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Model Parameters/Model_11';

N_MU = 200;
CT = zeros(1,5);
t2t = zeros(1,5);
twitch_amp = zeros(1,5);
p2p =zeros(52,5);
FR_half = zeros(N_MU,1);
Af = zeros(N_MU,52);
fusion = zeros(N_MU,52);
f_exp = [2 4 6 8 10 12.5 14 16 18 20 22 25 28 30 33.3 40 50 66.6 80 100];
force_exp = [2.5723 7.717 16.399 29.582 44.373 50.482 62.379 63.344 67.203 73.234 77.492 83.923 87.46 88.103 91.318 96.141 97.428 100.32 100 98.392];
[~,loc] = min(abs(force_exp-50));
f_half_exp = f_exp(loc);
fusion_exp = [-0.28169 17.465 36.62 58.31 74.93 85.915 90.141 92.394 94.93 96.056 96.62 97.465 97.465 97.746 98.028 97.746 97.465 97.465 97.746 96.62];


Lce_vec = [1 1.2 1.1 0.9 0.8];

for i = 100
    if i <= 147
        fiber_type = 'slow';
    else
        fiber_type = 'fast';
    end
    cd('/Users/akiranagamori/Documents/Github/Twitch-Based-Muscle-Model/Development_Code/Data')
    load(['Data_full_' num2str(i)])
    cd(code_folder)
    for j = 1:5
        Data = Data_all{j};
        CT(j) = Data{2,1};
        t2t(j) = Data{2,5};
        FR_half(i) = Data{2,6};
        
        twitch_amp(j) = Data{2,4};
        p2p(:,j) = Data{2,13};
        
        Af(i,:) = Data{2,10}';
        fusion(i,:) = Data{2,11}';
        
        f_eff = Data{2,9};
        Lce = Lce_vec(j);
        if strcmp(fiber_type,'slow')
            a_f = 0.56;
            n_f0 = 2.1;
            n_f1 = 5;
        elseif strcmp(fiber_type,'fast')
            a_f = 0.56;
            n_f0 = 2.1;
            n_f1 = 3.3;
        end
        n_f = n_f0 +n_f1* (1/Lce-1);
        Af_Song = 1-exp(-(f_eff./(a_f*n_f)).^n_f);
        
        figure(2)
        plot(Data{2,9},Data{2,10}*100,'LineWidth',1,'color',[11,19,43]/255)
        xlim([0 3])
        hold on
        plot(f_eff,Af_Song*100,'LineWidth',1,'color',[252,163,17]/255)
        
        figure(3)
        plot(Data{2,9},Data{2,11}*100,'LineWidth',1,'color',[11,19,43]/255)
        xlim([0 3])
        hold on
        
        figure(4)
        %ax_4 =
        plot(Data{2,10}*100,Data{2,11}*100,'LineWidth',1,'color',[11,19,43]/255)
        %ax_4.Color(4) = 0.5;
        hold on
        
        
    end
end

figure(2)
xlabel('Discharge Rate (f_{0.5})')
ylabel('Activation (%)')
set(gca,'TickDir','out');
set(gca,'box','off')

figure(3)
plot(f_exp./f_half_exp,fusion_exp,'LineWidth',1,'color',[252,163,17]/255)
xlabel('Discharge Rate (f_{0.5})')
ylabel('Fusion (%)')
set(gca,'TickDir','out');
set(gca,'box','off')
%legend('Lce = 1','Lce = 1.2', 'Lce = 1.1', 'Lce = 0.9', 'Lce = 0.8')

x_scale = 20/2.9;
y_scale =20/2.1;

x = [0.6 2.0 4.6 7 9.5 11.2 12.7 13.5 13.8];
y = [0.25 3 5 7.2 8.8 9.4 9.5 10.1 10.1];

Force = x*x_scale;
Fusion = y*y_scale;

figure(4)
hold on
plot(Force,Fusion,'LineWidth',2,'color',[252,163,17]/255)
xlabel('Activation (%)')
ylabel('Fusion (%)')
set(gca,'TickDir','out');
set(gca,'box','off')


figure(5)
plot(Data{2,9},1-p2p./twitch_amp(2),'LineWidth',1,'color',[11,19,43]/255)
xlim([0 3])
xlabel('Discharge Rate (f_{0.5})')
ylabel('Fusion (%)')
set(gca,'TickDir','out');
set(gca,'box','off')
