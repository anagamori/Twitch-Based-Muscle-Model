

Fs = 1000; %sampling frequency
time = 0:1/Fs:5; %simulation time
T_S = 0.043;
a_S1 = 1.76;
a_S2 = 0.96;
S = 0;
S_vec = zeros(1,length(time));
for t = 1:length(time)
    if Force(t) < 0.1
        a_S = a_S1;
    else
        a_S = a_S2;
    end
    S_dot = (a_S - S)/T_S;
    S = S_dot*1/Fs + S;
    
    S_vec(t) = S;
end

