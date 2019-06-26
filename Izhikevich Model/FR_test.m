function FR = FR_test(Fs,time,amp,parameter)

%input = [zeros(1,1*Fs) amp/2*[0:1/Fs:2] amp*ones(1,length(time)-1*Fs-length(amp*[0:1/Fs:2]))];
input = [zeros(1,1*Fs) amp*ones(1,length(time)-1*Fs)];
[~,binary] = Izhikevich(time,input,parameter,Fs);
spike_time = find(binary(3*Fs+1:end));
ISI = diff(spike_time)/(Fs/1000);
FR = mean(1./ISI*1000);

end