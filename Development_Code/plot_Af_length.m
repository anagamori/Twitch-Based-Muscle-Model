close all
clc
clear all

code_folder = '/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Development_Code';
data_folder = '/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Model Parameters/Model_11';

N_MU = 200;

f_exp = [2 4 6 8 10 12.5 14 16 18 20 22 25 28 30 33.3 40 50 66.6 80 100];
force_exp = [2.5723 7.717 16.399 29.582 44.373 50.482 62.379 63.344 67.203 73.234 77.492 83.923 87.46 88.103 91.318 96.141 97.428 100.32 100 98.392];
[~,loc] = min(abs(force_exp-50));
f_half_exp = f_exp(loc);
fusion_exp = [0 17.465 36.62 58.31 74.93 85.915 90.141 92.394 94.93 96.056 96.62 97.465 97.465 97.746 98.028 97.746 97.465 97.465 97.746 96.62];

Fs = 10000;

for i = 188
    
    cd('/Users/akiranagamori/Documents/Github/Twitch-Based-Muscle-Model/Development_Code/Data')
    load('modelParameter_v2');
    cd(code_folder)
    
    parameter = modelParameter.parameterMatrix(i,:);
    if i <= modelParameter.index_slow
        fiber_type = 'slow';
        FR_test_vec = [2:1:50 100];
    else
        fiber_type = 'fast';
        FR_test_vec = [2:1:90 100 120 150 200 300];
    end
    
    %
    Af = zeros(5,length(FR_test_vec));
    p2p = zeros(5,length(FR_test_vec));
    Lce = 1;
    tic
    [Data] = MUModel_test_full(parameter,Lce,0,FR_test_vec,fiber_type,0,0,Fs);
    Af(1,:) = Data{2,10};
    p2p(1,:) = Data{2,13};
    
    FR_half = Data{2,6};
    
    figure(1)
    plot(FR_test_vec,Data{2,10},'color',[11,19,43]/255)
    xlim([0 3])
    hold on
    
    figure(3)
    plot(Data{2,10}*100,Data{2,11}*100,'color',[11,19,43]/255)
    hold on
    
    f_eff = Data{2,9};
    
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
    
    figure(1)
    plot(FR_test_vec,Af_Song,'color',[252,163,17]/255)
    xlim([0 3])
    
    Lce_vec = [0.8 0.9 1.1 1.2];
    for j = 1:4
        Lce = Lce_vec(j);
        
        [Data] = MUModel_test_full(parameter,Lce,FR_half,FR_test_vec,fiber_type,0,0,Fs);
        
        Af(j+1,:) = Data{2,10};
        p2p(j+1,:) = Data{2,13};
        
        figure(1)
        plot(FR_test_vec,Data{2,10},'color',[11,19,43]/255)
        xlim([0 3])
        
        figure(3)
        plot(Data{2,10}*100,Data{2,11}*100,'color',[11,19,43]/255)
        hold on
        
        f_eff = Data{2,9};
        
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
        
        figure(1)
        plot(FR_test_vec,Af_Song,'color',[252,163,17]/255)
        xlim([0 3])
        
    end
    toc
end

figure(1)
set(gca,'TickDir','out');
set(gca,'box','off')

fusion = (1 - p2p./(p2p(5,1)))*100;
%%
figure(2)
plot(FR_test_vec,fusion,'color',[11,19,43]/255)
% hold on 
% plot(f_exp./f_half_exp,fusion_exp,'color',[252,163,17]/255)
xlim([0 3])
set(gca,'TickDir','out');
set(gca,'box','off')


figure(3)
plot(force_exp,fusion_exp,'color',[252,163,17]/255)
plot(0:1:100,0:1:100,'--k')
xlim([0 100])
set(gca,'TickDir','out');
set(gca,'box','off')

%%
% Lce = 1;
% FR_test_vec = [2 20 40 70 120];
% [Data] = MUModel_test_full(parameter,Lce,0,FR_test_vec,fiber_type,0,0,Fs);