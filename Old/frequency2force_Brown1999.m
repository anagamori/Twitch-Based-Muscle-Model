clc
close all
clc

Lce = 1;
Freq = [0.4555 0.9 1.1867 1.7742 3.5226];
Act = [0.1421, 0.4, 0.613, 0.818 1];

FR_half_reference = 14.4;
FR_reference = [3 5 8 10 15 20 30 40 50 60 80]/FR_half_reference;
Force_reference = [6.34 7.3 13.7 22.7 55.1 71.5 87.2, 93.7 96 98.9 100]/100;

figure(1)
plot(Freq,Act)
ylim([0 1])
hold on 
plot(FR_reference,Force_reference)