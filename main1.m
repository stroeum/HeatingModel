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

% flag.plot = 'on';
% U = 19.2;
% Heating1D(U,0,flag);
% Heating1D(U,2.5,flag);
% stop

flag.complex_ions = 1;
flag.complex_chem = 0;
flag.plot         = 'off';
flag.movie        = 'off';
U                 = [17.2 18.2 19.2 20.2 21.2 22.2 24.2];
% t_br_1p   = [];
% t_br_075p = [];
% 
%% Calculations
for ii=1:length(U)
    ii
    T0km(ii)  = 0;%Heating1D(U(ii),  0, flag); % 1-D
    T30km(ii) = 0;%Heating1D(U(ii), 30, flag); % 1-D
    T50km(ii) = 0;%Heating1D(U(ii), 50, flag); % 1-D
    T70km(ii) = 0;%Heating1D(U(ii), 70, flag); % 1-D
end

UL1p=[18.3 19.2 20.4 20.8 22 22.3 23.6 23.9 24.4]; % _kV
TL1p=[2.3 1.4 0.7 0.6 0.27 0.27 0.13 0.12 0.09]*1e-6; % _s
TC1p=10.^((UL1p-23.5)/(15.75-23.5)*-5+(UL1p-15.75)/(23.5-15.75)*-7); % _s

for ii=1:length(UL1p)
    ii
    t_br_1p(ii) = Heating1D( (UL1p(ii)-.2)/satm(0) + .2, 0, flag); % 1-D
end

UL075p=[17.9 18.8 19.2 20.8]; % _kV
TL075p=[0.22 0.18 0.13 0.05]*1e-6; % _s

for ii=1:length(UL075p)
    ii
    t_br_075p(ii) = Heating1D( (UL075p(ii)-.2)/satm(2.5) + .2, 2.5, flag); % 1-D
end

%% Plots
figure(1)
% semilogy(U, t_br_1p(:,1)*1e9,'r-', U, t_br_1p(:,2)*1e9,'r--', U, t_br_1p(:,3)*1e9,'r-.',UL1p,TL1p*1e9,'ks',UL1p,TC1p*1e9,'ko',U, t_br_075p(:,1)*1e9,'b-', U, t_br_075p(:,2)*1e9,'b--', U, t_br_075p(:,3)*1e9,'b-.',UL075p,TL075p*1e9,'ks')
semilogy(UL1p, t_br_1p*1e9,'r-',UL1p,TL1p*1e9,'ks',UL1p,TC1p*1e9,'ko',UL075p, t_br_075p*1e9,'b-',UL075p,TL075p*1e9,'ks')
legend('p = p_0','[Larsson et al., 1998]','[Cernak et al., 1995]','p = .75 p_0','Location','best');
legend('boxoff');
xlabel('U (kV)')
ylabel('\tau_{br} (ns)')
set(gca,'XMinorTick','on','YMinorTick','on','TickDir','out')
grid on

%%
figure(2)

subplot(221)
semilogy(U, 1e6*T0km, 'k-', U, 1e6*T30km, 'k-', U, 1e6*T50km, 'k-', U, 1e6*T70km, 'k-')
% axis([17 25 1e-6 2e3])
xlabel('E N_0/N (kV/cm)')
ylabel('Transition time \tau (\mus)')
set(gca,'XMinorTick','on','YMinorTick','on','TickDir','out')
grid on

subplot(222)
semilogy(U, 1e6*T0km, 'k-', U, 1e6*T30km*satm(30)^.5, 'k-', U, 1e6*T50km*satm(50)^.5, 'k-', U, 1e6*T70km*satm(70)^.5, 'k-')
% axis([17 25 1e-6 2e3])
xlabel('E N_0/N (kV/cm)')
ylabel('Transition time \tau (N/N_0)^.5 (\mus)')
set(gca,'XMinorTick','on','YMinorTick','on','TickDir','out')
grid on

subplot(223)
semilogy(U, 1e6*T0km, 'k-', U, 1e6*T30km*satm(30), 'k-', U, 1e6*T50km*satm(50), 'k-', U, 1e6*T70km*satm(70), 'k-')
axis([17 25 1e-4 2e1])
xlabel('E N_0/N (kV/cm)')
ylabel('Transition time \tau N/N_0 (\mus)')
set(gca,'XMinorTick','on','YMinorTick','on','TickDir','out')
grid off

subplot(224)
semilogy(U, 1e6*T0km, 'k-', U, 1e6*T30km*satm(30)^2, 'k-', U, 1e6*T50km*satm(50)^2, 'k-', U, 1e6*T70km*satm(70)^2, 'k-') 
axis([17 25 1e-4 2e1])
xlabel('E N_0/N (kV/cm)')
ylabel('Transition time \tau N^2/N_0^2 (\mus)')
set(gca,'XMinorTick','on','YMinorTick','on','TickDir','out')
grid off
