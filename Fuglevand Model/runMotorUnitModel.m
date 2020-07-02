close all
clear all
clc

Fs = 1000;
t = 0:1/Fs:15;
N_temp = 50:50:1000;
amp_temp = [0.025 0.05 0.1:0.1:1];
RP_temp = 10:10:150;


data_directory = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Data/Fuglevand/N_200_Ur_50';
code_directory = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Fuglevand Model';
%%
parpool(10)
maxForce = 3.3282e+04;
for k = 1:length(amp_temp)
    trialN = k; %+30; 
    % predefine model parameters
    amp = amp_temp(k); %0.15 for 0.05
    modelParameter.amp = amp;
    U = [zeros(1,1*Fs) (amp/2)*(0:1/Fs:2) amp*ones(1,length(t)-3*Fs-1)];
    
    modelParameter.N = 200;    
    modelParameter.RR = 17;     %65 for U_r = 0.8
    modelParameter.MFR = 8;   
    modelParameter.g_e = 1;    
    modelParameter.PFR1 = 35;   
    modelParameter.PFRD = 10;
    modelParameter.cv = 0.2;    
    modelParameter.RP = 100;    
    modelParameter.T_L = 90;    
    modelParameter.RT = 3;   
    
    Data = cell(1,10);
    tic
    parfor i = 1:10
        % Run motor unit model   
        tic
        output = MotorUnitModel_varCoV(t,U,modelParameter,Fs);       
        toc
        Data{i} = output;
    end
    toc
    
    cd (data_directory)
    save(['Trial_' num2str(trialN)],'Data','-v7.3')
    cd (code_directory)
    
    output_temp = Data{1};
    
    Force = output_temp.TotalForce(5*Fs+1:end)/maxForce;
    meanForce = mean(Force) 
    SD = std(Force)
    CoV = std(Force)/mean(Force)

    figure(1)
    plot(t,output_temp.TotalForce)
    xlabel('Time (s)','FontSize',14)
    ylabel('Force (AU)','FontSize',14)
    hold on
    
    k
    
    [pxx,f] = pwelch(Force-mean(Force),[],[],0:0.1:30,Fs,'power');
    figure(2)
    plot(f,pxx./sum(pxx)*100)
    xlabel('Frequency (Hz)','FontSize',14)
    ylabel('Power (%Total Power)','FontSize',14)
    hold on

end


