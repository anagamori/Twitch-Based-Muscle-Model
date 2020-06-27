%==========================================================================
% fit_MU_length.m
% Author: Akira Nagamori
% Last update: 6/23/2020
% Descriptions:
%   Randomly adjust model parameters to generate temporary sets of
%   parameters
%==========================================================================


close all
clc
clear all

code_folder = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Development_Code';
data_folder = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Model Parameters/Model_11';



simOpt = 0;
Fs = 5000;
parpool(10)
for m = 69:200
    m
    MU_No = m;
    if m <= 147
        MU_type = 'slow';
    else
        MU_type = 'fast';
    end
    cd('/Users/akira/Documents/Github/Twitch-Based-Muscle-Model/Development_Code/Data')
    load(['Data_v2_' num2str(MU_No)])
    cd(code_folder)
    
    parameter = Data{2,12};
    FR_half = Data{2,6};
    FR_temp = 0.5:0.1:3;
    FR_test = FR_temp*FR_half;
    
    Lce_vec = [1.2,1.1,0.9,0.8];
    
    phi_vec = 0.5:0.1:2;
    
    Data_all = cell(4,length(phi_vec));
    error_mat = zeros(4,length(phi_vec));
    
    %parpool(5)
    for j = 1:length(Lce_vec)
        tic
        Lce = Lce_vec(j);
        parfor i = 1:length(phi_vec)
            
            S = parameter(1);
            C = parameter(2);
            k_1 = parameter(3);
            k_2 = parameter(4);
            k_3 = parameter(5);
            k_4 = parameter(6);
            tau_1 = parameter(7);
            tau_2 = parameter(8);
            N = parameter(9);
            K = parameter(10);
            tau_3 = parameter(11);
            alpha = parameter(12);
            
            ptb = phi_vec(i);
            
            a_k_3 = k_3*ptb;
            a_N = -N*ptb;
            a_K = -K*ptb;
            a_alpha = alpha*ptb;
            
            k_3 = a_k_3*(Lce-1) + k_3; % a_k_3 > 0
            N = a_N*(Lce-1) + N; % a_k_3 < 0
            K = a_K*(Lce-1) + K;
            alpha = a_alpha*(Lce-1) + alpha; % a_k_3 > 0
            
            param = [S,C,k_1,k_2,k_3,k_4,tau_1,tau_2,N,K,tau_3,alpha];
            [Data_temp] = MUModel_test_length(param,Lce,FR_half,FR_test,MU_type,0,simOpt,Fs);
            hold on
            
            Data_all{j,i} = Data_temp;
            error_mat(j,i) = Data_temp{2,8};
        end
        toc
    end
    
    %%
    error_ecc = sum(error_mat(1:2,:));
    [~,loc_min_ecc] = min(error_ecc);
    
    error_con = sum(error_mat(3:4,:));
    [~,loc_min_con] = min(error_con);
    
    Data_length = cell(4,1);
    
    for n = 1:4
        if n <= 2
            temp = Data_all{n,loc_min_ecc};
        else
            temp = Data_all{n,loc_min_con};
        end
        
        f_eff = FR_temp;
        
        if strcmp(MU_type,'slow')
            a_f = 0.56;
            n_f0 = 2.1;
            n_f1 = 5;
        elseif strcmp(MU_type,'fast')
            a_f = 0.56;
            n_f0 = 2.1;
            n_f1 = 3.3;
        end
        n_f = n_f0 +n_f1* (1/Lce_vec(n)-1);
        Af_Song = 1-exp(-(f_eff./(a_f*n_f)).^n_f);
        
        figure(m)
        plot(temp{2,9},temp{2,10},'b','LineWidth',1)
        hold on
        plot(f_eff,Af_Song,'k')
        xlim([0 3])
        
        Data_length{n} = temp;
    end
    
    phi = [phi_vec(loc_min_ecc) phi_vec(loc_min_con)];
    cd('/Users/akira/Documents/Github/Twitch-Based-Muscle-Model/Development_Code/Data')
    save(['Data_length_v2_' num2str(MU_No)],'Data_length')
    save(['phi_v2_' num2str(MU_No)],'phi')
    cd(code_folder)
    
end