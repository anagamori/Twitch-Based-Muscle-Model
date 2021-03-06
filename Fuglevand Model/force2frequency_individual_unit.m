close all
clear all
clc

Fs = 1000;
time = 0:1/Fs:5;
time_spike = 0:1/Fs:3;

N = 200;
RR = 65;     %65 for U_r = 0.8
MFR = 8;
g_e = 1;
PFR1 = 35;
PFRD = 10;

i = 1:N; %motor unit identification index
a = log(RR)/N; %coefficient to establish a range of threshold values
RTE = exp(a*i); %recruitment threshold excitation
RTEn = exp(a*N); %recruitment threshold of the last motor unit
PFR = PFR1 - PFRD * (RTE./RTEn); %peak firing rate

RP = 100;
T_L = 90;
RT = 3;
b = log(RP)/N; %coefficient to establish a range of twich force values
P = exp(b*i); %force generated by a motor unit as a function of its recruitment threshold
c = log(100)/log(RT); %coefficient to establish a range of contraction time values
T = (T_L.* (1./P).^(1/c))./1000; %contraction time
t_twitch = 0:1/Fs:1;
twitch = zeros(N,length(t_twitch));

for j = 1:N
    twitch(j,:) =  P(j).*t_twitch./T(j).*exp(1-t_twitch./T(j));
end

f_test = 1:1:70;
cv = 0;
mean_force = zeros(1,length(f_test));
p2p_force = zeros(1,length(f_test));

mean_force_all = zeros(N,length(f_test));
fusion_all = zeros(N,length(f_test));
t2t = zeros(N,1);

for n = 1 %:200
    n
    for f = 1:length(f_test)
        spike_train = zeros(1,length(time));
        force = zeros(1,length(time));
        
        temp = spikeTrainGenerator(time_spike,Fs,f_test(f));
        spike_train(1*Fs+1:4*Fs+1) = temp;
        for t = 1:length(time)
            spike_train_temp = zeros(1,length(time));
            if spike_train(t) == 1
                if ~any(force)
                    spike_train_temp(t) = 1;
                    force_temp = conv(spike_train_temp,twitch(n,:));
                    force = force+ force_temp(1:length(time));
                else
                    spike_train_temp(t) = 1;
                    ISI = 1/f_test(f);
                    StimulusRate = T(n)/ISI;
                    if StimulusRate > 0 && StimulusRate <= 0.4
                        g = 1;
                    elseif StimulusRate > 0.4
                        S_MU = 1 - exp(-2*(StimulusRate)^3);
                        g = (S_MU/StimulusRate)/0.3;
                    end                    
                    force_temp = conv(spike_train_temp,g*twitch(n,:));
                    force = force+ force_temp(1:length(time));
                end
            end
        end
        mean_force(f) = mean(force(3.75*Fs:4*Fs));
        p2p_force(f) = max(force(3*Fs:4*Fs)) - min(force(3*Fs:4*Fs));
        fusion = 1-p2p_force./p2p_force(1);
        
    end
    
    figure(4)
    ax_4 = plot(mean_force./mean_force(end)*100,fusion*100,'color',[11,19,43]/255);
    ax_4.Color(4) = 0.5;
    hold on
    t2t(n) = p2p_force(1)/mean_force(end);
    mean_force_all(n,:) = mean_force;
    fusion_all(n,:) = fusion;
end

% save('mean_force_all','mean_force_all')
% save('fusion_all','fusion_all')
% save('t2t','t2t')