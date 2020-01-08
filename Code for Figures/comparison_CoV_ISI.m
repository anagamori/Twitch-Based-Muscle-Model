%==========================================================================
% comparison_CoV_ISI.m
% Author: Akira Nagamori
% Last update: 7/18/19
% Descriptions:
%   This code is used to generate fig XX in the paper
%==========================================================================
close all
clear all
clc

%%
code_folder = '/Users/akira/Documents/Github/Twitch-Based-Muscle-Model/Code for Figures';
figure_folder = '/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Figures';

Fs = 10000;
CoV = zeros(10,2);
for i = 1:3
    if i == 1
        condition = 'Model_8_no_CoV_50_Ur_Rec_3';
        
        data_folder = ['/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Data/withTendon/' condition];
        cd(data_folder)
        load('cov_Force')
        load('mean_pxx')
        cd(code_folder)
        
        color_code = [0 0 0];
    elseif i == 2
        condition = 'Model_8_10_CoV_50_Ur_Rec_3';
       
        data_folder = ['/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Data/withTendon/' condition];
        cd(data_folder)
        load('cov_Force')
        load('mean_pxx')
        cd(code_folder)
        
        color_code = [230 57 70]/255;
        
    elseif i == 3
        condition = 'Model_8_20_CoV_50_Ur_Rec_3';
       
        data_folder = ['/Users/akira/Documents/GitHub/Twitch-Based-Muscle-Model/Data/withTendon/' condition];
        cd(data_folder)
        load('cov_Force')
        load('mean_pxx')
        cd(code_folder)
        
        color_code = [37  65 178]/255;
    end
    
    %%
    CoV(:,i) = cov_Force(:,2);
    
    
end

figure(4)
%subplot(2,2,4)
boxplot(CoV,'Labels',{'CoV ISI = 0%','CoV ISI = 10%','CoV ISI = 20%'})
ylabel('CoV for Force (%)','FontSize',10)
%yticks(0.6:0.2:1.4)
set(gca,'TickDir','out');
set(gca,'box','off')
%ax = gca;
%ax.FontSize = 6;

% figure(4)
% fig = gcf;
% %linkaxes([ax1,ax2,ax3,ax4],'y')
% fig.PaperUnits = 'inches';
% fig.PaperPosition = [0 0 6.92 6.92];
