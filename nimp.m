%% Breakdown time as a function of the applied voltage
% U is applied between two electrodes separated by 1 cm

%% Update November 23, 2009
%
% # Creation of the main file
% # Plot of 1-D, 0-D cst p, 0-D cst n model results

%% Update November 24, 2009
%
% # Addition of the experimental results by Larsson [1998] to the plots

%% Initation
clear all
close all
clc

%% Parameters
flag.complex_ions = 1;
flag.complex_chem = 0;
flag.plot         = 'off';
flag.movie        = 'off';
U                 = [19.2 20.2 21.2 22.2 24.2];

%% Calculations
for ii=1:length(U)
    ii
    h=0;
    T0km(ii,1)  = Heating0Dn(U(ii), h, flag); % 0-D, n cst
    T0km(ii,2)  = Heating0Dp(U(ii), h, flag); % 0-D, p cst
  
    h=70;
    T70km(ii,1) = Heating0Dn(U(ii), h, flag); % 0-D, n cst
    T70km(ii,2) = Heating0Dp(U(ii), h, flag); % 0-D, p cst
end

%% Plots
figure(1)
semilogy(U, T0km(:,1)*1e6,'r--',U, T0km(:,2)*1e6,'g-.')
ylim([.1 1])
legend('0-D n cst','0-D p cst')
legend('boxoff');
xlabel('E N_0/N (kV/cm)')
ylabel('Transition time \tau (\mus)')
set(gca,'XMinorTick','on','YMinorTick','on','TickDir','out')
grid off

%%
figure(2)
semilogy(U, T70km(:,1)*1e3,'r--',U, T70km(:,2)*1e3,'g-.')
ylim([10 100])
legend('0-D n cst','0-D p cst')
legend('boxoff');
xlabel('E N_0/N (kV/cm)')
ylabel('Transition time \tau (ms)')
set(gca,'XMinorTick','on','YMinorTick','on','TickDir','out')