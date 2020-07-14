close all
clear all
clc

%%
code_folder = '/Users/akira/Documents/Github/Twitch-Based-Muscle-Model/Code for Figures';
figure_folder = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Figures';

amp_vec = 0.01:0.01:1;
for j = 2
    
    for i = 1:2
        if j == 1
            data_folder = '/Volumes/DATA2/PLOS_CB_Data/Fuglevand/singleUnit';
            if i == 1
                test_unit = 100;
            else
                test_unit = 180;
            end
            color_code = [100 100 100]/255;
        else
            
            if i == 1
                data_folder = '/Volumes/DATA2/PLOS_CB_Data/noTendon/Model_singleUnit';
                test_unit = 91;
            else
                data_folder = '/Volumes/DATA2/PLOS_CB_Data/noTendon/Model_singleUnit_constantCoV';
                test_unit = 91;
            end
            color_code = [230 57 70]/255;
        end
        cd(data_folder)
        load(['mean_Force_' num2str(test_unit)])
        load(['std_Force_' num2str(test_unit)])
        load(['std_Force_dt_' num2str(test_unit)])
        load(['cov_Force_' num2str(test_unit)])
        load(['cov_Force_dt_' num2str(test_unit)])
        load(['mean_FR_' num2str(test_unit)])
        load(['CoV_ISI_' num2str(test_unit)])
        cd(code_folder)
        
        mean_mean_Force = mean(mean_Force);
        mean_Force_norm = mean_Force./mean_mean_Force(end)*100;
        mean_Force_norm(mean_Force_norm==0) = nan;
        figure(1)
        shadedErrorBar(amp_vec*100,nanmean(mean_Force_norm),nanstd(mean_Force_norm),'lineprops',{'color',color_code,'LineWidth',2,'markerfacecolor',color_code})
        set(gca,'TickDir','out');
        set(gca,'box','off')
        hold on
        
        std_Force_norm = std_Force./mean_mean_Force(end)*100;
        std_Force_norm(std_Force_norm==0) = nan;
        
        figure(2)
        shadedErrorBar(amp_vec*100,nanmean(std_Force_norm),nanstd(std_Force_norm),'lineprops',{'color',color_code,'LineWidth',2,'markerfacecolor',color_code});
        set(gca,'TickDir','out');
        set(gca,'box','off')
        hold on
        
        
        figure(3)
        shadedErrorBar(amp_vec*100,nanmean(cov_Force),nanstd(cov_Force),'lineprops',{'color',color_code,'LineWidth',2,'markerfacecolor',color_code});
        set(gca,'TickDir','out');
        set(gca,'box','off')
        hold on
        
        mean_FR(mean_FR == 0) = nan;
        figure(4)
        shadedErrorBar(amp_vec*100,nanmean(mean_FR),nanstd(mean_FR),'lineprops',{'color',color_code,'LineWidth',2,'markerfacecolor',color_code});
        set(gca,'TickDir','out');
        set(gca,'box','off')
        hold on
        
        figure(5)
        shadedErrorBar(amp_vec*100,nanmean(CoV_ISI),nanstd(CoV_ISI),'lineprops',{'color',color_code,'LineWidth',2,'markerfacecolor',color_code});
        set(gca,'TickDir','out');
        set(gca,'box','off')
        hold on
        
    end
    
end

figure(1)
xlabel('Activation (%)','FontSize',14)
ylabel('Force (%Maximum)','FontSize',14)
ylim([0 105])
set(gca,'TickDir','out');
set(gca,'box','off')
legend('Fuglevand model','New model without tendon','New model with tendon','location','northwest')
% cd (figure_folder)
% saveas(gcf,'activation2meanForce_FV_comparison','pdf')
% cd (code_folder)

figure(2)
%xlabel('Mean Force (%)','FontSize',14)
%plot(x,x*a+b,'LineWidth',2,'Color','k')
xlabel('Mean Force (%Maximum Force)','FontSize',14)
ylabel('SD (%Maximum Force)','FontSize',14)
legend('Fuglevand model','New model without tendon','New model with tendon','location','northwest')
%yticks([0.05 0.1 0.15 0.2 0.25])
xlim([0 100])
set(gca,'TickDir','out');
set(gca,'box','off')
% cd (figure_folder)
% saveas(gcf,'activation2SD_FV_comparison','pdf')
% cd (code_folder)

figure(3)
xlabel('Mean Force (%Maximum Force)','FontSize',14)
ylabel('CoV (%)','FontSize',14)
xlim([0 100])
legend('Fuglevand model','New model without tendon','New model without tendon','location','northwest')
set(gca,'TickDir','out');
set(gca,'box','off')
% cd (figure_folder)
% saveas(gcf,'activation2CoV_FV_comparison','pdf')
% cd (code_folder)

figure(4)
xlabel('Synaptic Input (%Maximum)','FontSize',14)
ylabel('CoV of Force (%)','FontSize',14)
xlim([0 100])
set(gca,'TickDir','out');
set(gca,'box','off')

