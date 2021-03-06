% close all
% clear all
% clc

Fs = 10000; % sampling frequency of simulation
t = 0:1/Fs:5; % time vector

% simulated length change of muscle
% First column of Fig 3 in Mileusnic et al., (2006)
lengthChange = 0.01*sin(2*pi*10*t)+1; %[0.95.*ones(1,1*Fs) 0.11*t(1:1.1818*Fs)+0.95 1.08*ones(1,length(t)-1*Fs-1.1818*Fs)];
lengthChange = lengthChange + 0.01*sin(2*pi*2*t);
% smooth out the length change to avoid sharp transitions, which might
% affect differentiation 
Lce = lengthChange;  %smooth(lengthChange,100);
%Lce = Lce'; % muscle length
Vce = [0 diff(Lce).*Fs]; % muscle velocity
Ace = [0 diff(Vce).*Fs]; % muscle acceleration


% Fusimotor drive
gamma_dynamic = 20; % gamma dynamic
gamma_static = 20; % gamma static

% Run a simulation 
[outputIa,outputII,AP_bag1,AP_primary_bag2,AP_primary_chain] = muscleSpindle(Lce,Vce,Ace,Fs,gamma_dynamic,gamma_static);

%
U = outputIa/2000+outputII/2000;
Ia_delay = 15*Fs/1000;
II_delay = 15*Fs/1000;
T_U = 0.03;
U_eff = 0;
U_eff_2 = 0;
U_eff_vec = zeros(1,length(t));
U_eff_2_vec = zeros(1,length(t));
for i = 1:length(t) 
    if i > II_delay
    U_eff_dot = (outputIa(i-Ia_delay)/2000+outputII(i-II_delay)/2000-U_eff)/T_U;
    U_eff = U_eff_dot*1/Fs + U_eff;
    end
    U_eff_2_dot = (U_eff-U_eff_2)/T_U;
    U_eff_2 = U_eff_2_dot*1/Fs + U_eff_2;
    
    U_eff_vec(i) = U_eff;
    U_eff_2_vec(i) = U_eff_2;
end

%
[pxx,f] = pwelch(U_eff_2_vec-mean(U_eff_2_vec),gausswin(2*Fs),2*Fs*0.9,0:0.1:100,Fs,'power');
figure()
plot(f,pxx)

