close all
clear all
clc

code_folder = '/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Development_Code';
cd('/Users/akiranagamori/Documents/GitHub/Twitch-Based-Muscle-Model/Model Parameters/Model_11')
load('fusion')
load('Af')
load('FR_half')
cd(code_folder)

fusion_old = fusion;
Af_old = Af;

cd('/Users/akiranagamori/Documents/Github/Twitch-Based-Muscle-Model/Development_Code/Data')
load('fusion')
load('Af')
cd(code_folder)

fusion_new = fusion;
Af_new = Af;

x_scale = 20/2.9;
y_scale =20/2.1;

x = [0.6 2.0 4.6 7 9.5 11.2 12.7 13.5 13.8];
y = [0.25 3 5 7.2 8.8 9.4 9.5 10.1 10.1];

Force = x*x_scale;
Fusion = y*y_scale;

f_exp = [2 4 6 8 10 12.5 14 16 18 20 22 25 28 30 33.3 40 50 66.6 80 100];
force_exp = [2.5723 7.717 16.399 29.582 44.373 50.482 62.379 63.344 67.203 73.234 77.492 83.923 87.46 88.103 91.318 96.141 97.428 100.32 100 98.392];
[~,loc] = min(abs(force_exp-50));
f_half_exp = f_exp(loc);
fusion_exp = [0 17.465 36.62 58.31 74.93 85.915 90.141 92.394 94.93 96.056 96.62 97.465 97.465 97.746 98.028 97.746 97.465 97.465 97.746 96.62];

figure(1)
plot(mean(Af_old)*100,mean(fusion_old)*100,'LineWidth',1,'color',[11,19,43]/255)
hold on
plot(mean(Af_new)*100,mean(fusion_new)*100,'LineWidth',1,'color',[11,19,43]/255)
plot(Force,Fusion,'LineWidth',2,'color',[252,163,17]/255)
plot(force_exp,fusion_exp,'LineWidth',2,'color',[252,163,17]/255)
xlabel('Activation (%)')
ylabel('Fusion (%)')
set(gca,'TickDir','out');
set(gca,'box','off')