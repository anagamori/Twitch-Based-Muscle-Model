%close all
clear all
clc

N = 200;
RR = 110;     %65 for U_r = 0.8
MFR = 8;
g_e = 1.5;
PFR1 = 20;
PFRD = -30;

i = 1:N; %motor unit identification index
a = log(RR)/N; %coefficient to establish a range of threshold values
RTE = exp(a*i); %recruitment threshold excitation
RTEn = exp(a*N); %recruitment threshold of the last motor unit
PFR = PFR1 - PFRD * (RTE./RTEn); %peak firing rate
PFRn = PFR1 - PFRD; %peak firing rate of the last motor unit
Emax = RTEn + (PFRn - MFR)/g_e; %maximum excitatory input

E = 0:0.001:1;
E = E*Emax;

testUnit = 1;
FR = g_e.*(E - RTE(testUnit)) + MFR;
FR(FR>=PFR(testUnit)) = PFR(testUnit);
FR(FR<=MFR(testUnit)) = MFR(testUnit);

%%
figure(1)
plot(E/Emax,FR)
xlabel('Normalized Excitation')
ylabel('Discharge Rate (Hz)')
set(gca,'TickDir','out');
set(gca,'box','off')
hold on 

figure(2)
plot(E,FR)
xlabel('Excitation (AU)')
ylabel('Discharge Rate (Hz)')
set(gca,'TickDir','out');
set(gca,'box','off')
hold on 

figure(3)
plot(RTE/Emax,'o')
hold on 