%==========================================================================
% new_model_parameterFit_ST_length.m
% Author: Akira Nagamori
% Last update: 2/22/119
%==========================================================================

close all
clear all
clc

%% Folder name
code_folder = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model';
data_folder = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Data/FT';

%% 
for trialN = 231:300
param_Matrix_1 = zeros(10,11);
param_Matrix_2 = zeros(10,11);
param_Matrix_3 = zeros(10,11);
param_Matrix_4 = zeros(10,11);

cd(data_folder)
load(['Data_' num2str(trialN)])
param = Data{2,12};
cd(code_folder)
param_Matrix = Data{2,12};
k_3 = param_Matrix(:,5);
k_4 = param_Matrix(:,6);
N = param_Matrix(:,9);
K = param_Matrix(:,10);

for i = 1:10
    
    cd(data_folder)
    load(['Data_' num2str(trialN) '_' num2str(i+10)])
    cd(code_folder)    
    param_Matrix_1(i,:) = Data{2,12};
    error_1(i) = Data{2,8};
    cd(data_folder)
    load(['Data_' num2str(trialN) '_' num2str(i+20)])
    cd(code_folder)
    param_Matrix_2(i,:) = Data{2,12};
    error_2(i) = Data{2,8};
    cd(data_folder)
    load(['Data_' num2str(trialN) '_' num2str(i+30)])
    cd(code_folder)
    param_Matrix_3(i,:) = Data{2,12};
    error_3(i) = Data{2,8};
    cd(data_folder)
    load(['Data_' num2str(trialN) '_' num2str(i+40)])
    cd(code_folder)
    param_Matrix_4(i,:) = Data{2,12};
    error_4(i) = Data{2,8};
    
end


error_vec = [error_1 error_2 error_3 error_4];
index = find(error_vec>(quantile(error_vec,0.75)+3*iqr(error_vec)));
%% 
L0_long = [0.8 0.9 1.1 1.2];
L0_vec = [0.8*ones(10,1);0.9*ones(10,1);1.1*ones(10,1);1.2*ones(10,1)];
L0_vec(index) = [];
L0_vec = [L0_vec;1];

k_3_long = [param_Matrix_1(:,5) param_Matrix_2(:,5) param_Matrix_3(:,5) param_Matrix_4(:,5)]';
k_3_vec = reshape(k_3_long',[],1);
k_3_vec(index) = [];
k_3_vec = [k_3_vec;k_3];
p_k_3 = polyfit(L0_vec,k_3_vec,1)

k_4_long = [param_Matrix_1(:,6) param_Matrix_2(:,6) param_Matrix_3(:,6) param_Matrix_4(:,6)]';
k_4_vec = reshape(k_4_long',[],1);
k_4_vec(index) = [];
k_4_vec = [k_4_vec;k_4];
p_k_4 = polyfit(L0_vec,k_4_vec,1)

N_long = [param_Matrix_1(:,9) param_Matrix_2(:,9) param_Matrix_3(:,9) param_Matrix_4(:,9)]';
N_vec = reshape(N_long',[],1);
N_vec(index) = [];
N_vec = [N_vec;N];
p_N = polyfit(L0_vec,N_vec,1)

K_long = [param_Matrix_1(:,10) param_Matrix_2(:,10) param_Matrix_3(:,10) param_Matrix_4(:,10)]';
K_vec = reshape(K_long',[],1);
K_vec(index) = [];
K_vec = [K_vec;K];
p_K = polyfit(L0_vec,K_vec,1)

parameter = [param_Matrix(1:4) p_k_3(1) p_k_3(2) p_k_4(1) p_k_4(2) param_Matrix(7) param_Matrix(8) p_N(1) p_N(2) p_K(1) p_K(2) param_Matrix(11)];

cd(data_folder)
save(['MU_' num2str(trialN)],'parameter')
cd(code_folder)

end

% %%
% x = L0_vec;
% Y = [k_3_vec,k_4_vec,N_vec,K_vec];
% [n,d] = size(Y);
% X = cell(n,1);
% for i = 1:n
%     X{i} = [eye(d) x(i)*eye(d)];
% end
% 
% [beta,Sigma] = mvregress(X,Y,'algorithm','cwls')
