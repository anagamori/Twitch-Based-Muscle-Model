close all
clear all
clc

x_scale = 20/2.9;
y_scale =20/2.1;

x = [0.6 2.0 4.6 7 9.5 11.2 12.7 13.5 13.8];
y = [0.25 3 5 7.2 8.8 9.4 9.5 10.1 10.1];

Force = x*x_scale;
Fusion = y*y_scale;

figure
plot(Force,Fusion)
xlabel('Mean Force (%)','FontSize',14)
ylabel('Fusion (%)','FontSize',14)
set(gca,'TickDir','out');
set(gca,'box','off')