close all
clearvar
clc


cd('D:\MU_Population_Model\Old_model')
load('Data_12_1.mat')
cd

Force = output.ForceTendon(end-1*Fs+1:end);
mean(Force(end-1*Fs+1:end))/444.5844
CoV = std(Force)/mean(Force)
%        clear output
mean_FR = zeros(1,200);
CoV_FR = zeros(1,200);

for n = 1:200
    spike_time = find(output.spike_train(n,end-1*Fs+1:end));
    ISI = diff(spike_time)/(Fs/1000);
    mean_FR(n) = mean(1./ISI*1000);
    CoV_FR(n) = std(ISI)/mean(ISI)*100;
end

cd('D:\MU_Population_Model\SDN_test_IN_8')
load('Data_12_1.mat')