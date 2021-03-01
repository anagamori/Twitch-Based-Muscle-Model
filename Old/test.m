Fs = 10000;
time = 0:1/Fs:3;
tau_1 = 0.001;
tau_2 = 0.002;
R_temp_1 = 1-exp(-time./tau_1);
R_temp_2 = 1-exp(-time./tau_2);
x = R_temp_1.*R_temp_2;
y = zeros(1,length(time));
y(1*Fs) = 1;
tic 
z = conv(x,y,'full');
toc