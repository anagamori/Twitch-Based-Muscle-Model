
Fs = 2000; %sampling frequency
time = 0:1/Fs:5; %simulation time
f_env = zeros(1,length(time));
f_env(1*Fs:4*Fs) = 2;
%f_env(1*Fs) = 2;
Lce = 1;
% T_f_1 = 0.0343;
% T_f_2 = 0.0227;
% T_f_3 = 0.047;
% T_f_4 = 0.0252;
T_f_1 = 0.0206;
T_f_2 = 0.0136;
T_f_3 = 0.0282;
T_f_4 = 0.0151;
a_f = 0.56;
n_f_0 = 2.1;
n_f_1 = 3.3; %
n_f  = n_f_0+n_f_1*(1/Lce-1);

f_int_dot = 0;
f_int = 0;
f_eff_dot = 0;
f_eff = 0;

f_eff_vec = zeros(1,length(time));

for t = 1:length(time)
    Af = 1 - exp(-(f_eff/(a_f*n_f))^n_f );
    if f_eff_dot >= 0 
        T_f = T_f_1*Lce^2 + T_f_2*f_env(t);
    else
        T_f = (T_f_3 + T_f_4*Af)/Lce;
    end
    f_int_dot = (f_env(t) - f_int)/T_f;
    f_int = f_int_dot/Fs + f_int;
    f_eff_dot = (f_int - f_eff)/T_f;
    f_eff = f_eff_dot/Fs + f_eff;
    
    f_eff_vec(t) = f_eff;
end

figure(2)
plot(time,f_env/2,'k')
hold on 
plot(time,f_eff_vec/2)

