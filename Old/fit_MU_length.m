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

MU_type = 'slow';

simOpt = 0;
Fs = 5000;

for n = 6:10
    MU_No = n;
    cd('/Users/akira/Documents/Github/Twitch-Based-Muscle-Model/Development_Code/Data')
    load(['Data_' num2str(MU_No)])
    cd(code_folder)
    
    parameter = Data{2,12};
    FR_half = Data{2,6};
    FR_temp = 0.5:0.1:3;
    FR_test = FR_temp*FR_half;
    
    Lce_vec = [1.2,1.1,0.9,0.8];
    Data_all = cell(1,4);
    for j = 1:length(Lce_vec)
        j
        r = rand(4,20);
        r = r*0.5;
        error_vec = zeros(1,size(r,2));
        
        Lce = Lce_vec(j);
        
        for i = 1:size(r,2)
            i
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
            
            if Lce < 1
                k_3 = k_3 - k_3*r(1,i); %r(1);
                N = N + N*r(2,i); %r(3);
                alpha = alpha - alpha*r(3,i);
                K = K + K*r(4,i);
            else
                k_3 = k_3 + k_3*r(1,i); %r(1);
                N = N - N*r(2,i); %r(3);
                alpha = alpha + alpha*r(3,i);
                K = K - K*r(4,i);
            end
            param = [S,C,k_1,k_2,k_3,k_4,tau_1,tau_2,N,K,tau_3,alpha];
            
            
            [Data_temp] = MUModel_test_length(param,Lce,FR_half,FR_test,MU_type,0,simOpt,Fs);
            
            error_vec(i) = Data_temp{2,8};
            
            %         figure(2)
            %         plot(Data{2,9},Data{2,10},'LineWidth',1)
            %         xlim([0 3])
            %         hold on
        end
        
        %%
        %     if strcmp(MU_type,'slow')
        %         a_f = 0.56;
        %         n_f0 = 2.1;
        %         n_f1 = 5;
        %     elseif strcmp(MU_type,'fast')
        %         a_f = 0.56;
        %         n_f0 = 2.1;
        %         n_f1 = 3.3;
        %     end
        %     n_f = n_f0 +n_f1* (1/Lce-1);
        %     Af_Song = 1-exp(-(FR_temp ./(a_f*n_f)).^n_f);
        %     figure(2)
        %     plot(FR_temp,Af_Song,'k','LineWidth',2)
        
        %%
        [val,idx] = min(error_vec);
        
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
        
        if Lce < 1
            k_3 = k_3 - k_3*r(1,idx); %r(1);
            N = N + N*r(2,idx); %r(3);
            alpha = alpha - alpha*r(3,idx);
            K = K + K*r(4,idx);
        else
            k_3 = k_3 + k_3*r(1,idx); %r(1);
            N = N - N*r(2,idx); %r(3);
            alpha = alpha + alpha*r(3,idx);
            K = K - K*r(4,idx);
        end
        
        k_3
        N
        alpha
        K
        
        param = [S,C,k_1,k_2,k_3,k_4,tau_1,tau_2,N,K,tau_3,alpha];
        [Data] = MUModel_test_length(param,Lce,FR_half,FR_test,MU_type,1,simOpt,Fs);
        hold on
        
        Data_all{j} = Data;
    end
    
    cd('/Users/akira/Documents/Github/Twitch-Based-Muscle-Model/Development_Code/Data')
    save(['Data_length_' num2str(MU_No)],'Data_all')
    cd(code_folder)
    
end