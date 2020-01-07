
Fs = 1000;
time = 0:1/Fs:10;

T = 0.03;
t_twitch = 0:1/Fs:1;
twitch = t_twitch./T.*exp(1-t_twitch./T);
twitch = twitch./(sum(twitch));

input =(rand(1,length(time))-0.5)*2;
output = conv(input,twitch);
output = output(1:length(input));
figure(1)
plot(time,input)
hold on
plot(time,output)

%%
[pxx_input,f] = pwelch(input,[],[],0:0.1:100,Fs);
[pxx_output,~] = pwelch(output,[],[],0:0.1:100,Fs);
gain = pxx_output./pxx_input;

figure(2)
semilogx(f,db(gain))

%%
H = tf([exp(1)/T],[1,2/T,1/(T)^2]);
figure(3)
h = bodeplot(H);
p = getoptions(h);
p.FreqUnits = 'Hz';
p.MagUnits = 'db';
setoptions(h,p);