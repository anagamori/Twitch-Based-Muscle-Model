close all
clear all
clc

data_directory = '/Volumes/DATA2/PLOS_CB_Data/Fuglevand/singleUnit';
code_directory = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Fuglevand Model/';
%load ('Input')
Fs = 1000;
t = 0:1/Fs:10;
nTrial = 10;
amp_vec = 0.01:0.01:1;

pxxAll = zeros(length(amp_vec),301);
Force_mat = zeros(10,15*Fs+1);

mean_Force = zeros(nTrial,length(amp_vec));
std_Force = zeros(nTrial,length(amp_vec));
cov_Force = zeros(nTrial,length(amp_vec));
std_Force_dt = zeros(nTrial,length(amp_vec));
cov_Force_dt = zeros(nTrial,length(amp_vec));
mean_FR = zeros(nTrial,length(amp_vec));
CoV_ISI = zeros(nTrial,length(amp_vec));
mean_pxx = zeros(length(amp_vec),301);

test_unit = 180;

for j = 1:length(amp_vec)
    trialN = j; 
    % predefine model parameters
   
    pxx = zeros(10,301);
    cd (data_directory)
    load(['Trial_' num2str(test_unit) '_' num2str(trialN+1)])
    cd (code_directory)
    for i = 1:10
        output = Data{i};
        Force = output.TotalForce(5*Fs+1:end);
        window = 1*Fs;
        bp = [1:window:length(Force)];
        Force_dt = detrend(Force,1,bp);
        
        mean_Force(i,j) = mean(Force);
        std_Force(i,j) = std(Force);
        std_Force_dt(i,j) = std(Force_dt);
        cov_Force(i,j) =  std_Force(i,j)/mean_Force(i,j)*100;
        cov_Force_dt(i,j) =  std(Force_dt)/mean_Force(i,j)*100;
        
        [pxx(i,:),f] = pwelch(Force-mean(Force),gausswin(5*Fs),2.5*Fs,0:0.1:30,Fs,'power');
        pxx(i,:) = pxx(i,:)./sum(pxx(i,:));
        pxx(i,:) = smooth(pxx(i,:),10);
        
         spike_time = find(output.SpikeTrain(5*Fs+1:end));
        ISI = diff(spike_time)/(Fs/1000);
        mean_FR(i,j) = mean(1./ISI*1000);
        CoV_ISI(i,j) = std(ISI)/mean(ISI)*100;
        
      
    end
    
    mean_pxx(j,:) = mean(pxx);
    
        figure(11)
    plot(Force)
    hold on
    
%     cd(data_directory)
%     save(['Force_mat_' num2str(j)],'Force_mat')
%     cd(code_directory)
    
end

cd(data_directory)
save(['mean_Force_' num2str(test_unit)],'mean_Force')
save(['std_Force_' num2str(test_unit)],'std_Force')
save(['std_Force_dt_' num2str(test_unit)],'std_Force_dt')
save(['cov_Force_' num2str(test_unit)],'cov_Force')
save(['cov_Force_dt_' num2str(test_unit)],'cov_Force_dt')
save(['mean_FR_' num2str(test_unit)],'mean_FR')
save(['CoV_ISI_' num2str(test_unit)],'CoV_ISI')
cd(code_directory)


maxForce = mean(mean_Force(:,end));

meanForce_plot = mean(mean_Force)./mean(mean_Force(:,end));
figure(1)
plot([0 amp_vec*100],[0 meanForce_plot*100])
xlabel('Excitatory Drive (%)')
ylabel('Force (%)')

figure(2)
plot(amp_vec*100,mean(cov_Force))
xlabel('Excitatory Drive (%)')
ylabel('CoV (%)')

figure(3)
plot(amp_vec*100,mean(std_Force./maxForce*100))
xlabel('Excitatory Drive (%)')
ylabel('SD (AU)')

figure(4)
plot(f,mean_pxx)

% cd (dataFolder)
% save('CoV_300','CoVAll')
% save('pxx_300','pxxAll')
% cd (codeFolder)