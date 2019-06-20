%==========================================================================
% build_MNs.m
% Author: Akira Nagamori
% Last update: 6/18/19
%==========================================================================

close all
clear all
clc

code_folder = '/Users/akira/Documents/Github/Twitch-Based-Muscle-Model/Izhikevich Model';
model_parameter_folder =  '/Users/akira/Documents/Github/Twitch-Based-Muscle-Model/Model Parameters/Model_4_Ur_50_constantT2T';

%%
Fs = 10000;
time = 0:1/Fs:5;

%%
for n = 86
    cd(model_parameter_folder )
    load('modelParameter')
    cd(code_folder)
    MDR = modelParameter.FR_half/2;
    PDR = modelParameter.FR_half*2;
    U_th = modelParameter.U_th_new;
    testUnit = n;
    
    U_th_target = U_th(testUnit)*100;
    MDR_target = MDR(testUnit);
    PDR_target = PDR(testUnit);
    
    %% Fixed parameters
    parameter.a = 0.02;
    parameter.c = -65;
    parameter.d = 100;
    parameter.v = -65;
    parameter.alpha = 0.04;
    parameter.beta = 5;
    parameter.gamma = 140;
    
    %%
    for m = 1:3
        for i = 1:3
            if i == 1
                b_vec = -0.2:0.001:0.6;
                FR_mat = zeros(length(b_vec),2);
                for j = 1:length(b_vec)
                    %% Model parameters
                    
                    % Parameter to be searched
                    parameter.b = b_vec(j);
                    
                    %% Input
                    for k = 1:2
                   
                        if k == 1
                            amp = U_th_target*0.01;
                            input = [zeros(1,1*Fs) amp/2*[0:1/Fs:2] amp*ones(1,length(time)-1*Fs-length(amp*[0:1/Fs:2]))];
                        else
                            amp = U_th_target;
                            input = [zeros(1,1*Fs) amp/2*[0:1/Fs:2] amp*ones(1,length(time)-1*Fs-length(amp*[0:1/Fs:2]))];
                        end
                        %% Run Izhikevich model
                        [v_vec,binary] = Izhikevich(time,input,parameter,Fs);
                        spike_time = find(binary(3*Fs+1:end));
                        ISI = diff(spike_time)/(Fs/1000);
                        FR_mat(j,k) = mean(1./ISI*1000);
                    end
                end
                b_index_1 = isnan(FR_mat(:,1));
                b_index_2 = ~isnan(FR_mat(:,2));
                b_index = find(b_index_1==1&b_index_2==1,1,'first');
                parameter.b = b_vec(b_index);
            elseif i == 2
                a_vec = 0:0.001:0.1;
                FR_mat = zeros(length(a_vec),1);
                for j = 1:length(a_vec)
                    %% Model parameters
                    parameter.a = a_vec(j);
                    %% Input
                    amp = U_th_target;
                    input = [zeros(1,1*Fs) amp/2*[0:1/Fs:2] amp*ones(1,length(time)-1*Fs-length(amp*[0:1/Fs:2]))];
                    
                    %% Run Izhikevich model
                    [v_vec,binary] = Izhikevich(time,input,parameter,Fs);
                    spike_time = find(binary(3*Fs+1:end));
                    ISI = diff(spike_time)/(Fs/1000);
                    FR_mat(j) = mean(1./ISI*1000);
                end
                [error,a_index] = min(abs(MDR_target-FR_mat));
                parameter.a = a_vec(a_index);
                
            elseif i == 3
                d_vec = 1:300;
                FR_mat = zeros(length(d_vec),2);
                for j = 1:length(d_vec)
                    %% Model parameters
                    parameter.d = d_vec(j);
                    %% Input
                    for k = 1:2
                        if k == 1
                            amp = U_th_target;
                            input = [zeros(1,1*Fs) amp/2*[0:1/Fs:2] amp*ones(1,length(time)-1*Fs-length(amp*[0:1/Fs:2]))];
                        elseif k== 2
                            amp = 100;
                            input = [zeros(1,1*Fs) amp/2*[0:1/Fs:2] amp*ones(1,length(time)-1*Fs-length(amp*[0:1/Fs:2]))];
                        end
                        %% Run Izhikevich model
                        [v_vec,binary] = Izhikevich(time,input,parameter,Fs);
                        spike_time = find(binary(3*Fs+1:end));
                        ISI = diff(spike_time)/(Fs/1000);
                        FR_mat(j,k) = mean(1./ISI*1000);
                        
                    end
                    
                end
                [error,d_index] = min(abs(PDR_target-FR_mat(:,2))+ abs(MDR_target-FR_mat(:,1)));
                parameter.d = d_vec(d_index)
            end
        end
        
        input_vec = 0.5:0.1:100; %[0.5:0.1:0.9 1:100];
        for i = 1:length(input_vec)
            %% Input
            amp = input_vec(i);
            input = [zeros(1,1*Fs) amp/2*[0:1/Fs:2] amp*ones(1,length(time)-1*Fs-length(amp*[0:1/Fs:2]))];


            %% Run Izhikevich model
            [v_vec,binary] = Izhikevich(time,input,parameter,Fs);
            
            spike_time = find(binary(3*Fs+1:end));
            ISI = diff(spike_time)/(Fs/1000);
            mean_FR(i) = mean(1./ISI*1000);
            
        end
        
        [error_I,I_index] = min(abs(PDR_target-mean_FR));
        parameter.I_max = input_vec(I_index);
        parameter.I_min = U_th_target;
        
        figure(1)
        plot(input_vec,mean_FR)
        hold on
    end
    cd([model_parameter_folder '/MN'])
    save(['MN_' num2str(n)],'parameter')
    cd(code_folder)
end