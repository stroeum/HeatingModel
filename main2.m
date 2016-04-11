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
U                 = [17.2 18.2 19.2 20.2 21.2 22.2 24.2];

%% Calculations
for ii=1:length(U)
    ii
    h=0;
    flag.complex_ions = 1;
    T0km(ii,3)  = Heating1D( U(ii), h, flag); % 1-D
    
    flag.complex_ions = 0;
    T0km(ii,1)  = Heating0Dn(U(ii), h, flag); % 0-D, n cst
    T0km(ii,2)  = Heating0Dp(U(ii), h, flag); % 0-D, p cst
    T0km(ii,4)  = Heating1D( U(ii), h, flag); % 1-D
    
    h=70;
    flag.complex_ions = 1;
    T70km(ii,3) = Heating1D( U(ii), h, flag); % 1-D
    flag.complex_ions = 0;
    T70km(ii,1) = Heating0Dn(U(ii), h, flag); % 0-D, n cst
    T70km(ii,2) = Heating0Dp(U(ii), h, flag); % 0-D, p cst
    T70km(ii,4) = Heating1D( U(ii), h, flag); % 1-D
end

%% Plots
figure(1)
semilogy(U, T0km(:,1)*1e6,'r--',U, T0km(:,2)*1e6,'g-.',U, T0km(:,3)*1e6,'b-',U, T0km(:,4)*1e6,'k-')
% ylim([.1 10])
legend('0-D n cst','0-D p cst', '1-D w/ C ions', '1-D w/o C ions');
legend('boxoff');
xlabel('E N_0/N (kV/cm)')
ylabel('Transition time \tau (\mus)')
set(gca,'XMinorTick','on','YMinorTick','on','TickDir','out')
grid off

%%
figure(2)
semilogy(U, T70km(:,1)*1e3,'r--',U, T70km(:,2)*1e3,'g-.',U, T70km(:,3)*1e3,'b-',U, T70km(:,4)*1e3,'k-')
% ylim([10 1000])
legend('0-D n cst','0-D p cst', '1-D w/ C ions', '1-D w/o C ions');
legend('boxoff');
xlabel('E N_0/N (kV/cm)')
ylabel('Transition time \tau (ms)')
set(gca,'XMinorTick','on','YMinorTick','on','TickDir','out')
grid off
grid off

dTdtT=nut*sigma*E^2*2/(5*n*1e6*kB); %deg/s
dTdtVT=QVT*2/5*11600; %deg/s
dTdtVTO=QVTO*2/5*11600; %deg/s

loglog(t, dTdtT, '-', t, dTdtVT, '--', t, dTdtVTO, '-')
axis([1e-9/nno tf 1e7*nno*nno 1e12*nno*nno])
xlabel('Time (sec)', 'FontSize', 20)
ylabel('dT/dt (deg/s)', 'FontSize', 20)
legend( 'Fast heating rate', 'VT heating rate', 'VT heating rate due to quenchig by O', 0);
