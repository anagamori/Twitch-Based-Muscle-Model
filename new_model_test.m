
function [Data] = new_model_test(parameter,Lce,FR_half_temp,fiber_type)
%==========================================================================
% model_test.m
% Author: Akira Nagamori
% Last update: 2/25/19
% Descriptions
%   Run twitch and sweep simulations with a given set of parameters and
%   generate associated data with figures
%   This code REQUIRES parameters related to length dependence on Af
%   relationship (e.g., parameter(5), parameter(6), etc)
%==========================================================================
%% Simulation parameters
Fs = 1000; %sampling frequency
time = 0:1/Fs:5; %simulation time

S = parameter(1); %7;
C = parameter(2); %1.025;
k_1 = parameter(3); %14.625;
k_2 = parameter(4); %4.9375;
k_3 = parameter(5)*Lce + parameter(6); %17.41*Lce - 2.85;
k_4 = parameter(7)*Lce + parameter(8); %-7.67*Lce + 14.92;
tau_1 = parameter(9); %0.0051;
tau_2 = parameter(10); % 0.04;
N = parameter(11)*Lce + parameter(12); %-2.26*Lce + 4.20;
K = parameter(13)*Lce + parameter(14); %-0.044*Lce + 0.080;
alpha = parameter(15); %4.475;

%% 
for i = 1:2
    if i == 1
        % Generate a pulse
        FR_test = 1;
    elseif i == 2
        % Generate a set of spike trains at multiple frequencies
        FR_test = [2 5 8 10 12 15 18 20 25 30 40 50 60 70 80 100 200]; %10:10:100];
    end
    %==========================================================================
    % initialization
    mean_exc = zeros(1,length(FR_test));
    p2p_exc = zeros(1,length(FR_test));
    for f = 1:length(FR_test)
        
        %% Generate spike train
        FR = FR_test(f);
        spike = zeros(1,length(time));
        
        if i == 1
            % Generate a pulse
            spike(1*Fs) = 1;
        else
            % Generate spike train
            temp =spikeTrainGenerator(0:1/Fs:3,Fs,FR);
            spike(1*Fs:4*Fs) = temp;
        end
        
        %%
        % Model parameter initilization
        x = 0; % free calcium concentration
        y = 0; % concentraction of calcium bound to troponin
        z = 0; % porportion of bidning sites that formed cross-bridges (0-1)
        act = 0;
        
        x_vec = zeros(1,length(time));
        y_vec = zeros(1,length(time));
        z_vec = zeros(1,length(time));
        act_vec = zeros(1,length(time));
        R_temp = exp(-time/tau_1);
        R = zeros(1,length(time));
        
        %% Run a simulation
        for t = 1:length(time)
            %% Stage 1
            % Calcium diffusion to sarcoplasm
            spike_temp = zeros(1,length(time));
            if spike(t) == 1
                spike_temp(t) = 1;
                %temp = conv(spike_temp,R_temp);
                % R to depend on the normalized firing rate
                temp = conv(spike_temp,R_temp*(1+2*act^alpha));
                %             R_i = 1 + (5-1)*exp(-(1/FR)/0.1);
                %             temp = conv(spike_temp,R_temp*R_i);
                R = R + temp(1:length(time));
            end
            
            %% Stage 2
            % Calcium kinetics
            x_dot = k_1*(C-x-y)*R(t) - k_2*x*((C*(S-1))+x+y)-(k_3*x+k_4*y)*(1-y);
            y_dot = (1-y)*(k_3*x-k_4*y);
            x = x_dot/Fs + x;
            y = y_dot/Fs + y;
            
            if y < 0
                y_temp = 0;
            else
                y_temp = y;
            end
            y_int = y_temp^N/(y_temp^N+K^N);
            
            %% Stage 3
            % Dynamics of cross-bridge formation
            z_dot = (y_int-z)/tau_2;
            z = z_dot/Fs + z;
            act = z; % Muscle activation
            
            %% Store variables
            %x_vec(t) = x(t);
            x_vec(t) = x;
            y_vec(t) = y;
            z_vec(t) = z;
            act_vec(t) = act;
            
        end
        
        mean_exc(f) = mean(act_vec(3.5*Fs:4*Fs));
        p2p_exc(f) = max(act_vec(3.5*Fs:4*Fs))-min(act_vec(3.5*Fs:4*Fs));
        
        if i == 1
            %------------------------------------------------------------------
            % Twitch analysis
            %------------------------------------------------------------------
            % Find the peak amplitude of twitch and its location
            [pks,locs_0_100] = max(act_vec);
            % Display peak amplitude value
            pks
            % Time it takes from zero force to the peak
            t_0_100 = (locs_0_100-1*Fs)*1000/Fs % convert it in the unit of ms
            
            % Find a half peak amplitude
            peak_half = pks/2;
            % Find time from peak to half maximum
            [~,t_100_50] = min(abs(act_vec(locs_0_100:end)-peak_half));
            
            % Find the location of t_100_50
            locs_100_50 = locs_0_100+t_100_50;
            % Time it takes from peak to half maximum
            t_100_50 = t_100_50*1000/Fs
            
            % Find value of 40% of peak amplitude
            peak_40 = pks*0.4;
            % Find value of 10% of peak amplitude
            peak_10 = pks*0.1;
            
            % Find time when force reaches 40% maximum after peak
            [~,t_40] = min(abs(act_vec(locs_0_100:end)-peak_40));
            % Find time when force reaches 10% maximum after peak
            [~,t_10] = min(abs(act_vec(locs_0_100:end)-peak_10));
            
            % Time it takes from 40% peak to 10% peak
            t_40_10 = (t_10-t_40)*1000/Fs
            
            % Locatino of t_40_10
            locs_40_10 = locs_0_100+t_10;
            
            T = t_0_100/Fs;
            P =  pks/T;
            twitch_Milner_temp = P.*time.*exp(1-time/T);
            twitch_Milner = conv(spike,twitch_Milner_temp);
        end
        
        figure(i)
        plot(time,act_vec,'LineWidth',1)
        hold on
        if i == 1
            plot(time,twitch_Milner(1:length(time)))
            hold on
            plot([time(locs_0_100) time(locs_0_100)],[0 1],'--')
            text(time(locs_0_100),pks,num2str(t_0_100))
            hold on
            plot([time(locs_100_50) time(locs_100_50)],[0 1],'--')
            plot([time(locs_40_10) time(locs_40_10)],[0 1],'--')
            text(time(locs_100_50),peak_half,num2str(t_100_50))
            text(time(locs_40_10),peak_10,num2str(t_40_10))
            text(time(locs_0_100)-time(locs_0_100)*0.4,pks,num2str(pks))
        end
        
    end
    
    figure(i)
    xlabel('Time (s)','FontSize',14)
    ylabel('Activation','FontSize',14)
            
    if i == 2
        %% Sweep simulation  
        twitch2tetanus_ratio = p2p_exc(1)/mean_exc(f)
        fusion = 1-p2p_exc/p2p_exc(1);
        
        FR_new = 0.1:0.1:FR_test(end);
        Af_new = spline(FR_test,mean_exc,FR_new);
        
        if FR_half_temp == 0
            [~,loc] = min(abs(Af_new-0.5));
            FR_half = FR_new(loc);
        else
            FR_half = FR_half_temp
        end
        f_eff = FR_new/FR_half;
        
        if fiber_type == 'slow'
            a_f = 0.56;
            n_f0 = 2.1;
            n_f1 = 5;
        elseif fiber_type == 'fast'
            a_f = 0.56;
            n_f0 = 2.1;
            n_f1 = 3.3;
        end
        n_f = n_f0 +n_f1* (1/Lce-1);
        Af_Song = 1-exp(-(f_eff./(a_f*n_f)).^n_f);        
        
        %% Calculate error between the desired and generated activation-frequency relationship
        for j = 1:length(find(f_eff<3.5))
            error_temp(j) = abs(Af_Song(j)-Af_new(j));
        end
        error = sum(error_temp)
        
        figure(6)
        plot(FR_test/FR_half,mean_exc,'LineWidth',2,'color','b')
        xlabel('Frequency (f_{0.5})','FontSize',14)
        ylabel('Activation','FontSize',14)
        set(gca,'TickDir','out');
        set(gca,'box','off')
        hold on
        %plot(f_eff,Af_new,'color','r')
        plot(f_eff,Af_Song,'color','k','LineWidth',1)
        xlim([0 3])
        legend('New','Song')
        
        figure(7)
        plot(FR_test/FR_half,fusion,'LineWidth',2,'color','b')
        xlabel('Frequency (f_{0.5})','FontSize',14)
        ylabel('Fusion','FontSize',14)
        set(gca,'TickDir','out');
        set(gca,'box','off')
        hold on
        xlim([0 3])
        
        figure(8)
        plot(mean_exc./max(mean_exc),fusion,'LineWidth',2,'color','b')
        xlabel('Activation','FontSize',14)
        ylabel('Fusion','FontSize',14)
        set(gca,'TickDir','out');
        set(gca,'box','off')
        hold on
        plot(0:0.1:1,0:0.1:1,'--','color','k')
    end
    
end

Data = cell(2,12);
Data(1,:) = {'t_0_100','t_100_50','t_40_10','P_t','t2t','FR_half','Error_CT','Error_Af','Test frequencies','Af','Fusion','Parameters'};
Data{2,1} = t_0_100;
Data{2,2} = t_100_50;
Data{2,3} = t_40_10;
Data{2,4} = pks;
Data{2,5} = twitch2tetanus_ratio;
Data{2,6} = FR_half;
Data{2,7} = 0;
Data{2,8} = error;
Data{2,9} = FR_test/FR_half;
Data{2,10} = mean_exc;
Data{2,11} = fusion;
Data{2,12} = [S,C,k_1,k_2,k_3,k_4,tau_1,tau_2,N,K,alpha];

    function spikeTrain = spikeTrainGenerator(t,Fs,freq)
        
        spikeTrain = zeros(1,length(t));
        ISI = round(1/freq*Fs);
        numSpikes = round(length(t)/ISI);
        index = [1:numSpikes]*ISI;
        index(index>length(t)) = [];
        spikeTrain(index) = 1;
        spikeTrain(1) = 1;
        
    end

    function [var1,var2,var3,var4,var5,var6,var7,var8,var9,var10,var11] = parameter_Assigning(x)
        var1 = x(1);
        var2 = x(2);
        var3 = x(3);
        var4 = x(4);
        var5 = x(5);
        var6 = x(6);
        var7 = x(7);
        var8 = x(8);
        var9 = x(9);
        var10 = x(10);
        var11 = x(11);
    end

end