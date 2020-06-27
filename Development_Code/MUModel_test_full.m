
function [Data] = MUModel_test_full(parameter,Lce,FR_half_temp,FR_test_vec,fiber_type,plotOpt,simOpt,Fs)
%==========================================================================
% MUModel_test_full.m
% Author: Akira Nagamori
% Last update: 6/25/2020
% Descriptions:
%==========================================================================
%% Simulation parameters
time = 0:1/Fs:6; %simulation time

S = parameter(1); 
C = parameter(2); 
k_1 = parameter(3); 
k_2 = parameter(4);
k_3 = parameter(5);
k_4 = parameter(6);
k_4_i = k_4;
tau_1 = parameter(7); 
tau_2 = parameter(8);
N = parameter(9); 
K = parameter(10);
tau_3 = parameter(11);
alpha = parameter(12);
phi_1 = parameter(13);
phi_2 = parameter(14);

if Lce >= 1 
    k_3 = (phi_1*k_3)*(Lce-1) + k_3; % a_k_3 > 0
    N = (-phi_1*N)*(Lce-1) + N; % a_k_3 < 0
    K = (-phi_1*K)*(Lce-1) + K;
    alpha = (phi_1*alpha)*(Lce-1) + alpha; % a_k_3 > 0
elseif Lce < 1
    k_3 = (phi_2*k_3)*(Lce-1) + k_3; % a_k_3 > 0
    N = (-phi_2*N)*(Lce-1) + N; % a_k_3 < 0
    K = (-phi_2*K)*(Lce-1) + K;
    alpha = (phi_2*alpha)*(Lce-1) + alpha; % a_k_3 > 0
end

%%
if simOpt == 0
    nTrial = 2;
else
    nTrial = 1;
end
%%
for i = 1:nTrial
    if i == 1
        % Generate a pulse
        FR_test = 1;
    elseif i == 2
        FR_test = FR_test_vec;  
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
            temp =spikeTrainGenerator(0:1/Fs:4,Fs,FR);
            spike(1*Fs:5*Fs) = temp;
        end
        
        %%  initialization
        c = 0; % free calcium concentration
        cf = 0; % concentraction of calcium bound to troponin
        A = 0; % muscle activation
        
        c_vec = zeros(1,length(time));
        cf_vec = zeros(1,length(time));
        A_tilda_vec = zeros(1,length(time));
        A_vec = zeros(1,length(time));
        
        
        R_temp = 1-exp(-time/tau_1);
        R_temp_2 = exp(-time/tau_2);
        R = zeros(1,length(time));
        %%
        for t = 1:length(time)
            %
            % Cooperativity
            k_4 = k_4_i/(1+alpha*A);
            %% Stage 1
            % Calcium diffusion to sarcoplasm
            spike_temp = zeros(1,length(time));
            if spike(t) == 1
                spike_temp(t) = 1;
                temp = conv(spike_temp,R_temp_2.*R_temp);
                R = R + temp(1:length(time));
            end
       
            %%
            c_dot = k_1*(C-c-cf)*R(t) - k_2*c*(S-C+c+cf) - (k_3*c-k_4*cf)*(1-cf);
            cf_dot = (1-cf)*(k_3*c-k_4*cf);
            c = c_dot/Fs + c;
            cf = cf_dot/Fs + cf;
       
            %% Stage 2
            % Cooperativity and saturation
            if cf < 0
                cf_temp = 0;
            else
                cf_temp = cf;
            end
            A_tilda = cf_temp^N/(cf_temp^N+K^N);
            
            %% Stage 3
            % First-order dynamics to muscle activation, A
            A_dot = (A_tilda-A)/tau_3;
            A = A_dot/Fs + A;
            
            %% Store variables
            c_vec(t) = c;
            cf_vec(t) = cf;
   
            A_tilda_vec(t) = A_tilda;
            A_vec(t) = A;
            
            
        end
        
        mean_exc(f) = mean(A_vec(4.5*Fs:5*Fs));
        p2p_exc(f) = max(A_vec(4.5*Fs:5*Fs))-min(A_vec(4.5*Fs:5*Fs));
        
        if i == 1
            %------------------------------------------------------------------
            % Twitch analysis
            %------------------------------------------------------------------
            % Find the peak amplitude of twitch and its location
            [pks,locs_0_100] = max(A_vec);
            % Display peak amplitude value
            pks;
            % Time it takes from zero force to the peak
            t_0_100 = (locs_0_100-1*Fs)*1000/Fs; % convert it in the unit of ms
            
            % Find a half peak amplitude
            peak_half = pks/2;
            % Find time from peak to half maximum
            [~,t_100_50] = min(abs(A_vec(locs_0_100:end)-peak_half));
            
            % Find the location of t_100_50
            locs_100_50 = locs_0_100+t_100_50;
            % Time it takes from peak to half maximum
            t_100_50 = t_100_50*1000/Fs;
            
            % Find value of 40% of peak amplitude
            peak_40 = pks*0.4;
            % Find value of 10% of peak amplitude
            peak_10 = pks*0.1;
            
            % Find time when force reaches 40% maximum after peak
            [~,t_40] = min(abs(A_vec(locs_0_100:end)-peak_40));
            % Find time when force reaches 10% maximum after peak
            [~,t_10] = min(abs(A_vec(locs_0_100:end)-peak_10));
            
            % Time it takes from 40% peak to 10% peak
            t_40_10 = (t_10-t_40)*1000/Fs;
            
            % Locatino of t_40_10
            locs_40_10 = locs_0_100+t_10;
            
            T = t_0_100/1000;
            P =  pks/T;
            twitch_Milner_temp = P.*time.*exp(1-time/T);
            twitch_Milner = conv(spike,twitch_Milner_temp);
        end        
    end
    
    
    
    if i == 2
        %% Sweep simulation
        twitch2tetanus_ratio = p2p_exc(1)/mean_exc(f);
        fusion = 1-p2p_exc/p2p_exc(1);
        
        FR_new = 1:0.1:FR_test(end);
        Af_new = makima(FR_test,mean_exc,FR_new);
        
        if FR_half_temp == 0
            [~,loc] = min(abs(Af_new-0.5));
            FR_half = FR_new(loc);
        else
            FR_half = FR_half_temp;
        end
        f_eff = FR_test_vec/FR_half;
        
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
        error = sum(abs(Af_Song-mean_exc));
%         for j = 1:length(FR_test)
%             
%         end
%         error = sum(error_temp);
        if plotOpt == 1
            figure(6)
            plot(FR_test/FR_half,mean_exc,'LineWidth',2)
            %plot(FR_test,mean_exc,'LineWidth',2,'color','b')
            %xlim([0 100])
            xlabel('Frequency (f_{0.5})','FontSize',14)
            ylabel('Activation','FontSize',14)
            set(gca,'TickDir','out');
            set(gca,'box','off')
            hold on
            %plot(FR_new/FR_half,Af_new,'color','r')
            plot(f_eff,Af_Song,'color','k','LineWidth',1)
            xlim([0 3])
            legend('New','Song')
            movegui('northeast')
            
            figure(7)
            plot(FR_test/FR_half,fusion,'LineWidth',2)
            xlabel('Frequency (f_{0.5})','FontSize',14)
            ylabel('Fusion','FontSize',14)
            set(gca,'TickDir','out');
            set(gca,'box','off')
            hold on
            xlim([0 3])
            movegui('southwest')
            
            figure(8)
            plot(mean_exc,fusion,'LineWidth',2)
            xlabel('Activation','FontSize',14)
            ylabel('Fusion','FontSize',14)
            set(gca,'TickDir','out');
            set(gca,'box','off')
            hold on
            plot(0:0.1:1,0:0.1:1,'--','color','k')
            movegui('southeast')
        end
    end
    
end
if nTrial == 2
    Data = cell(2,13);
    Data(1,:) = {'t_0_100','t_100_50','t_40_10','P_t','t2t','FR_half','Error_CT','Error_Af','Test frequencies','Af','Fusion','Parameters','p2p'};
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
    Data{2,12} = [S,C,k_1,k_2,k_3,k_4,tau_1,tau_2,N,K,tau_3,alpha,phi_1,phi_2];
    Data{2,13} =  p2p_exc;
else
    Data = cell(2,13);
    Data(1,:) = {'t_0_100','t_100_50','t_40_10','P_t','t2t','FR_half','Error_CT','Error_Af','Test frequencies','Af','Fusion','Parameters','p2p'};
    Data{2,1} = t_0_100;
    Data{2,2} = t_100_50;
    Data{2,3} = t_40_10;
    Data{2,4} = pks;
end
    function spikeTrain = spikeTrainGenerator(t,Fs,freq)
        
        spikeTrain = zeros(1,length(t));
        ISI = round(1/freq*Fs);
        numSpikes = round(length(t)/ISI);
        index = [1:numSpikes]*ISI;
        index(index>length(t)) = [];
        spikeTrain(index) = 1;
        spikeTrain(1) = 1;
        
    end

end