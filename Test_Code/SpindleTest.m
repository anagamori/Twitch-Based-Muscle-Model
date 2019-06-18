%close all
%clear all
clc
%
% dataFolder = '/Users/akira/Documents/Github/Afferented-Muscle/CommonDrive/CommonDriveData';
% cd (dataFolder)
% load('output_Ia6000_II10000_Ib1000_GammaDynamic25_GammaStatic25_PIC_FB_1')
% cd ..
% cd ..
% Lce = output.Lce;
% Vce = output.Vce;
% Ace = output.Ace;

Fs = 10000;
Fs_Feedback = 10000;
%t = 0:1/Fs:(length(Lce)-1)/Fs;

t = 0:1/Fs:5;
lengthChange = [0.9.*ones(1,1*Fs) 0.18*t(1:1*Fs)+0.9 1.08.*ones(1,length(t)-1*Fs-1*Fs)];
%lengthChange = 0.0015*cos(1*2*pi.*t)+0.995;
% Lce = smooth(lengthChange,100);
%Lce = lengthChange;
%Lce = Lce';
Lce = output.Lce;
Vce = output.Vce;
Ace = output.Ace;

Lce = Lce(1);
Vce = 0;
Ace = 0;

[b,a] = butter(4,50/(Fs/2),'low');
[b100,a100] = butter(4,100/(Fs/2),'low');

for j = 1
    
    gamma_dynamic = 10;
    gamma_static = 10;
    
    figure(1)
    subplot(2,1,1)
    plot(t,Lce)
    subplot(2,1,2)
    plot(t,Vce)
    
    %% Feedback system parameters
    %% Spindle Model
    p = 2;
    R = 0.46; %length dependency of the force-velocity relationship
    a = 0.3;
    K_SR = 10.4649;
    K_PR = 0.15;
    M = 0.0002; % intrafusal fiber mass
    
    LN_SR = 0.0423;
    LN_PR = 0.89;
    
    L0_SR = 0.04; %in units of L0
    L0_PR = 0.76; %polar region rest length
    
    L_secondary = 0.04;
    X = 0.7;
    
    tau_bag1 = 0.149;
    freq_bag1 = 60;
    
    beta0_bag1 = 0.0605;
    beta1 = 0.2592;
    Gamma1 = 0.0289;
    
    G_bag1 = 20000;  %7000
    
    tau_bag2 = 0.205;
    freq_bag2 = 60;
    
    beta0 = 0.0822;
    beta2 = -0.046;
    Gamma2 = 0.0636;
    
    G_bag2 = 10000; %7250 %3800
    
    freq_chain = 90;
    
    beta0 = 0.0822;
    beta2_chain = - 0.069;
    Gamma2_chain = 0.0954;
    
    
    
    G_chain = 10000; %7250    %3000
    
    f_dynamic = 0;
    f_static = 0;
    T_ddot_bag1 = 0;
    T_dot_bag1 = 0;
    T_bag1 = 0;
    T_ddot_bag2 = 0;
    T_dot_bag2 = 0;
    T_bag2 = 0;
    T_ddot_chain = 0;
    T_dot_chain = 0;
    T_chain = 0;
    
    
    count = 1;
    
    for i = 1 %:length(Lce)
        L = Lce(i);
        L_dot = Vce(i);
        L_ddot = Ace(i);
        
        
        if Vce(i) >= 0
            C = 1;
        else
            C = 0.42;
        end
        
        df_dynamic = (gamma_dynamic^p/(gamma_dynamic^p+freq_bag1^p)-f_dynamic)/tau_bag1;
        f_dynamic = 1/Fs_Feedback*df_dynamic + f_dynamic;
        
        beta_bag1 = beta0_bag1 + beta1 * f_dynamic;
        Gamma_bag1 = Gamma1 * f_dynamic;
        T_ddot_bag1 = K_SR/M * (C * beta_bag1 * sign(L_dot-T_dot_bag1/K_SR)*((abs(L_dot-T_dot_bag1/K_SR))^a)...
            *(L-L0_SR-T_bag1/K_SR-R)+K_PR*(L-L0_SR-T_bag1/K_SR-L0_PR)+M*L_ddot+Gamma_bag1-T_bag1);
        temp_2(i) = T_ddot_bag1;
        %     if i > 5
        %         temp(i) = K_SR/M * (C * beta_bag1 * sign(L_dot-T_dot_bag1/K_SR)*((abs(L_dot-T_dot_bag1/K_SR))^a)...
        %                             *(L-L0_SR-T_bag1/K_SR-R)+K_PR*(L-L0_SR-T_bag1/K_SR-L0_PR)+M*L_ddot+Gamma_bag1-T_bag1);
        %         temp_filt(i) = (b100(5)*temp(i-4) + b100(4)*temp(i-3) + b100(3)*temp(i-2) + b100(2)*temp(i-1) + b100(1)*temp(i) ...
        %             - a100(5)*temp_filt(i-4) - a100(4)*temp_filt(i-3) - a100(3)*temp_filt(i-2) - a100(2)*temp_filt(i-1))/a100(1);
        %         T_ddot_bag1 = temp_filt(i);
        %     else
        %         temp(i) = K_SR/M * (C * beta_bag1 * sign(L_dot-T_dot_bag1/K_SR)*((abs(L_dot-T_dot_bag1/K_SR))^a)...
        %                             *(L-L0_SR-T_bag1/K_SR-R)+K_PR*(L-L0_SR-T_bag1/K_SR-L0_PR)+M*L_ddot+Gamma_bag1-T_bag1);
        %         temp_filt(i) = temp(i);
        %         T_ddot_bag1 = temp_filt(i);
        %     end
        T_dot_bag1 = T_ddot_bag1*1/Fs_Feedback + T_dot_bag1;
        T_bag1 = T_dot_bag1*1/Fs_Feedback + T_bag1;
        
        element1(i) = beta_bag1;
        element2(i) = sign(L_dot-T_dot_bag1/K_SR);
        element3(i) = (abs(L_dot-T_dot_bag1/K_SR))^a;
        element4(i) = L-L0_SR-T_bag1/K_SR-R;
        element5(i) = C;
        element6(i) = K_PR*(L-L0_SR-T_bag1/K_SR-L0_PR);
        element7(i) = M*L_ddot;
        element8(i) = Gamma_bag1;
        element9(i) = T_bag1;
        
        element10(i) = T_bag1;
        
        AP_bag1 = G_bag1*(T_bag1/K_SR-(LN_SR-L0_SR));
        
        df_static = (gamma_static^p/(gamma_static^p+freq_bag2^p)-f_static)/tau_bag2;
        f_static = 1/Fs_Feedback*df_static + f_static;
        
        beta_bag2 = beta0 + beta2 * f_static;
        Gamma_bag2 = Gamma2 * f_static;
        
        T_ddot_bag2 = K_SR/M * (C * beta_bag2 * sign(L_dot-T_dot_bag2/K_SR)*((abs(L_dot-T_dot_bag2/K_SR))^a)...
            *(L-L0_SR-T_bag2/K_SR-R)+K_PR*(L-L0_SR-T_bag2/K_SR-L0_PR)+M*L_ddot+Gamma_bag2-T_bag2);
        T_dot_bag2 = T_ddot_bag2*1/Fs_Feedback + T_dot_bag2;
        T_bag2 = T_dot_bag2*1/Fs_Feedback + T_bag2;
        
        AP_primary_bag2 = G_bag2*(T_bag2/K_SR-(LN_SR-L0_SR));
        AP_secondary_bag2 = G_bag2*(X*L_secondary/L0_SR*(T_bag2/K_SR-(LN_SR-L0_SR))+(1-X)*L_secondary/L0_PR*(L-T_bag2/K_SR-L0_SR-LN_PR));
        
        
        f_static_chain = gamma_static^p/(gamma_static^p+freq_chain^p);
        
        beta_chain = beta0 + beta2_chain * f_static_chain;
        Gamma_chain = Gamma2_chain * f_static;
        
        T_ddot_chain = K_SR/M * (C * beta_chain * sign(L_dot-T_dot_chain/K_SR)*((abs(L_dot-T_dot_chain/K_SR))^a)...
            *(L-L0_SR-T_chain/K_SR-R)+K_PR*(L-L0_SR-T_chain/K_SR-L0_PR)+M*L_ddot+Gamma_chain-T_chain);
        T_dot_chain = T_ddot_chain*1/Fs_Feedback + T_dot_chain;
        T_chain = T_dot_chain*1/Fs_Feedback + T_chain;
        
        AP_primary_chain = G_chain*(T_chain/K_SR-(LN_SR-L0_SR));
        AP_secondary_chain = G_chain*(X*L_secondary/L0_SR*(T_chain/K_SR-(LN_SR-L0_SR))+(1-X)*L_secondary/L0_PR*(L-T_chain/K_SR-L0_SR-LN_PR));
        
        
        S = 0.156; % 0.156
        
        if AP_bag1 < 0
            AP_bag1 = 0;
        end
        
        if AP_primary_bag2 < 0
            AP_primary_bag2 = 0;
        end
        
        if AP_primary_chain < 0
            AP_primary_chain = 0;
        end
        
        
        if AP_secondary_bag2 < 0
            AP_secondary_bag2 = 0;
        end
        
        if AP_secondary_chain < 0
            AP_secondary_chain = 0;
        end
        
        
        if AP_bag1 > (AP_primary_bag2+AP_primary_chain)
            Larger = AP_bag1;
            Smaller = AP_primary_bag2+AP_primary_chain;
        elseif AP_bag1 < (AP_primary_bag2+AP_primary_chain)
            Larger = AP_primary_bag2+AP_primary_chain;
            Smaller = AP_bag1;
        elseif AP_bag1 == (AP_primary_bag2+AP_primary_chain)
            Larger = 0;
            Smaller = 0;
        end
        OutputPrimary = Larger + S * Smaller;
        OutputSecondary = AP_secondary_bag2 + AP_secondary_chain;
        
        if OutputPrimary < 0
            OutputPrimary = 0;
        elseif OutputPrimary > 100000
            OutputPrimary = 100000;
        end
        if OutputSecondary < 0
            OutputSecondary = 0;
        elseif OutputSecondary > 100000
            OutputSecondary = 100000;
        end
        outputIa(i) = OutputPrimary;
        outputII(i) = OutputSecondary;
        %     if i == count
        %         outputIa(i) = OutputPrimary;
        %         outputII(i) = OutputSecondary;
        %         count = 10+count;
        %     else
        %         outputIa(i) = outputIa(i-1);
        %         outputII(i) = outputII(i-1);
        %     end
        actionPotential_bag1(i) = AP_bag1;
        actionPotential_bag2(i) = AP_primary_bag2;
        actionPotential_chain(i) = AP_primary_chain;
    end
    
    
    % figure(1)
    % subplot(3,1,1)
    % plot(t,outputIa)
    % subplot(3,1,2)
    % plot(t,outputII)
    % subplot(3,1,3)
    % plot(t,Lce)
    % hold on
    
    figure(2)
    subplot(4,1,1)
    plot(t,Lce,'b')
    xlabel('Time (s)')
    ylabel('Lce (L0)')
    subplot(4,1,2)
    plot(t,Vce,'b')
    xlabel('Time (s)')
    ylabel('Vce (L0/s)')
    
    subplot(4,1,3)
    plot(t,actionPotential_bag1)
    xlabel('Time (s)')
    ylabel('AP_{bag1}')
    
    hold on
    subplot(4,1,4)
    plot(t,actionPotential_bag2+actionPotential_chain)
    xlabel('Time (s)')
    ylabel('AP_{bag2} + AP_{chain}')
    
    hold on
    
    % temp = K_SR/M.*(element1.*element2.*element3.*element4.*element5+element6+element7+element8-element9);
    % temp_filt = filtfilt(b,a,temp);
    % figure(3)
    % ax(1) =  subplot(2,1,1);
    % plot(t,temp_filt,'Parent',ax(1));
    % hold on
    % ax(2) = subplot(2,1,2);
    % plot(t,outputIa,'Parent',ax(2))
    % hold on
    % linkaxes(ax,'x');
    
    % element1(i) = beta_bag1;
    %     element2(i) = sign(L_dot-T_dot_bag1/K_SR);
    %     element3(i) = (abs(L_dot-T_dot_bag1/K_SR))^a;
    %     element4(i) = L-L0_SR-T_bag1/K_SR-R;
    %     element5(i) = C;
    %     element6(i) = K_PR*(L-L0_SR-T_bag1/K_SR-L0_PR);
    %     element7(i) = M*L_ddot;
    %     element8(i) = Gamma_bag1;
    %     element9(i) = T_bag1;
    %
    %     element10(i) = T_bag1;
    
    % figure(4)
    % ax(1) = subplot(3,3,1);
    % plot(t,element1,'Parent',ax(1))
    % ylabel('\beta');
    % xlim([1.5 2.5])
    % hold on
    % ax(2) = subplot(3,3,2);
    % plot(t,element2,'Parent',ax(2))
    % hl = ylabel('$$sign(\dot{L}-\dot{T}/K^{SR})$$');
    % set(hl, 'Interpreter', 'latex');
    % xlim([1.5 2.5])
    % hold on
    % ax(3) = subplot(3,3,3);
    % plot(t,filtfilt(b,a,element3),'Parent',ax(3))
    % hl =  ylabel('$abs(\dot{L}-\dot{T}/K^{SR})^{a}$');
    % set(hl, 'Interpreter', 'latex');
    % xlim([1.5 2.5])
    % hold on
    % ax(4) = subplot(3,3,4);
    % plot(t,element4,'Parent',ax(4))
    % hl = ylabel('$L-L^{SR}_{0}-T/K^{SR}-R$');
    % set(hl, 'Interpreter', 'latex');
    % xlim([1.5 2.5])
    % hold on
    % ax(5) = subplot(3,3,5);
    % plot(t,element5,'Parent',ax(5))
    % hl = ylabel('C');
    % set(hl, 'Interpreter', 'latex');
    % xlim([1.5 2.5])
    % hold on
    % ax(6) = subplot(3,3,6);
    % plot(t,element1.*element2.*element3.*element4.*element5,'Parent',ax(6))
    % hl = ylabel({'$$C*\beta*sign(\dot{L}-\dot{T}/K^{SR})$$','$$\ *abs(\dot{L}-\dot{T}/K^{SR})^{a}$$','$$\ *(L-L^{SR}_{0}-T/K^{SR}-R)$$'});
    % set(hl, 'Interpreter', 'latex');
    % xlim([1.5 2.5])
    % hold on
    % ax(7) = subplot(3,3,7);
    % plot(t,element6+element7+element8-element9,'Parent',ax(7))
    % hl =  ylabel({'$$K^{SR}*(L-L^{SR}_{0}$$','$$\ -T/K^{SR}-L^{SR}_{0})+$$','$$\ M*\ddot{L}+\Gamma-T$$'});
    % set(hl, 'Interpreter', 'latex');
    % xlim([1.5 2.5])
    % hold on
    % ax(8) = subplot(3,3,8);
    % plot(t,filtfilt(b,a,element1.*element2.*element3.*element4.*element5+element6+element7+element8-element9),'Parent',ax(8))
    % ylabel('Eqaution (6)')
    % xlim([1.5 2.5])
    % hold on
    % ax(9) = subplot(3,3,9);
    % plot(t,outputIa,'Parent',ax(9))
    % ylabel('Ia Afferent Firing')
    % xlim([1.5 2.5])
    % hold on
    % linkaxes(ax,'x');
    
    % subplot(2,1,2)
    % plot(t,trapz(K_SR/M * (element1.*element1.*element3.*element4*C+element5+element6+element7+element8-element9)))
    % hold on
    
    % figure(4)
    % plot(t,element10)
    % hold on
    %
    % figure(5)
    % plot(t,element3)
    % hold on
    
    
end
%legend('Gamma Dynamic = 50 pps','Gamma Dynamic = 100 pps')

%hold on
%plot(t,output.Ia,'r')



