%close all
clc
clear all

code_folder = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Development_Code';
data_folder = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Model Parameters/Model_11';


simOpt = 0;
Fs = 10000;
%parpool(4)

for n = 1:200
    MU_No = n;
    if n <= 147
        MU_type = 'slow';
    else
        MU_type = 'fast';
    end
    cd('/Users/akira/Documents/Github/Twitch-Based-Muscle-Model/Development_Code/Data')
    load(['Data_' num2str(MU_No)])
    load(['phi_' num2str(MU_No)])
    cd(code_folder)
    
    parameter = [Data{2,12} phi];
    
    Lce_vec = [1.2,1.1,0.9,0.8];
    
    parameter(12) = 0.8;
    parameter(11) = parameter(11)+0.002;
    
    FR_test_vec = [2:2:100 200 300];
    %FR_test_vec = [2 4 6 8 10 15 20 25 30 50 100 200 300];
    
    Data_all = cell(1,5);
    
    Lce = 1;
    tic
    [Data] = MUModel_test_full(parameter,Lce,0,FR_test_vec,MU_type,1,simOpt,Fs);
    Data_all{1} = Data;
    FR_half = Data{2,6};
    
    parfor j = 1:length(Lce_vec)
        Lce = Lce_vec(j);
        [Data_temp] = MUModel_test_full(parameter,Lce,FR_half,FR_test_vec,MU_type,1,simOpt,Fs);
        Data_all{j+1} = Data_temp;
    end
    toc
    
    for i = 1:4
        
        temp = Data_all{i+1};
        f_eff = FR_test_vec./FR_half;
        
        if strcmp(MU_type,'slow')
            a_f = 0.56;
            n_f0 = 2.1;
            n_f1 = 5;
        elseif strcmp(MU_type,'fast')
            a_f = 0.56;
            n_f0 = 2.1;
            n_f1 = 3.3;
        end
        n_f = n_f0 +n_f1* (1/Lce_vec(i)-1);
        Af_Song = 1-exp(-(f_eff./(a_f*n_f)).^n_f);
        
        figure(6)
        plot(temp{2,9},temp{2,10},'b','LineWidth',1)
        hold on
        plot(f_eff,Af_Song,'k')
        xlim([0 3])
        
    end
    
    cd('/Users/akira/Documents/Github/Twitch-Based-Muscle-Model/Development_Code/Data')
    save(['Data_v2_' num2str(MU_No)],'Data')
    cd(code_folder)
end