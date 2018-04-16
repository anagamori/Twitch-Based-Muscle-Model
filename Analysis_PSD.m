close all
clear all
clc

dataFolder = '/Users/akiranagamori/Documents/GitHub/SDN Data';
codeFolder = '/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model';
amp_vec = 0.1:0.1:1;
test_amp = [1,3,5,8];
Fs = 40000;
Fs_2 = 1000;

mean_pxx_all_1 = zeros(length(test_amp),1001);
mean_pxx_all_2 = zeros(length(test_amp),1001);

freq = 0:0.1:100;

for k = 1:length(test_amp)
    i = test_amp(k);
    cd (dataFolder)
    load(['Force_FDI_' num2str(i)],'Force')
    Force_1 = Force;
    load(['Force_FDI_noTendon_' num2str(i)],'Force')
    Force_2 = Force;
    cd (codeFolder)
    
    pxx_Force_1 = zeros(10,1001); 
    pxx_Force_2 = zeros(10,1001); 
    
    
    for j = 1:10
        Force_temp = Force_1(j,:);        
        [pxx_Force_1(j,:),~] = pwelch(Force_temp-mean(Force_temp),[],[],0:0.1:100,Fs);
        %pxx_Force_1(j,:) = pxx_Force_1(j,:)./sum(pxx_Force_1(j,:));
        Force_temp_2 = Force_2(j,:);        
        [pxx_Force_2(j,:),~] = pwelch(Force_temp_2-mean(Force_temp_2),[],[],0:0.1:100,Fs_2);
        %pxx_Force_2(j,:) = pxx_Force_2(j,:)./sum(pxx_Force_2(j,:));
        
    end
    
    mean_pxx_all_1(k,:) = mean(pxx_Force_1);
    mean_pxx_all_2(k,:) = mean(pxx_Force_2);
   
    
end

%%
figure(1)
plot(freq,mean_pxx_all_1,'LineWidth',2);
xlim([0 60])
ylim([0 0.005])
yticks([0:0.001:0.05])
xlabel('Frequency (Hz)','FontSize',14,'fontweight','bold')
ylabel('Power (N^2)','FontSize',14,'fontweight','bold')
legend('10% MVC','30% MVC','50% MVC','80% MVC')

%%
figure(2)
plot(freq,mean_pxx_all_2,'LineWidth',2)
xlim([0 60])
ylim([0 0.005])
yticks([0:0.001:0.005])
xlabel('Frequency (Hz)','FontSize',14,'fontweight','bold')
ylabel('Power (N^2)','FontSize',14,'fontweight','bold')
legend('10% MVC','30% MVC','50% MVC','80% MVC')