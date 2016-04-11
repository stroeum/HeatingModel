%% August 28, 2003 Chemical heating of streamers
% This is a code to simulate heating in the streamer channel including effects of electron detachment by metastable O2
% References:
%
% * Lowke, J. Phys. D: Appl. Phys., 25, 202, 1992
% * Naidis, J. Phys. D: Appl. Phys., 32, 2649, 1999
% * Vidal et al., IEEE Trans. Plasma Sci, 30, 1339, 2002
% * Pasko, Theoretical modeling of sprites and jets, in Sprites, Elves and
% Intense Lightning Discharges (NATO Science Series II: Mathematics,
% Physics and Chemistry Vol. 225), ed. M. Fullekrug et al., Heidleberg:
% Springer, p. 253-293, 2006. for scaling.

%% Update April 7, 2007
%
% # We implement nested functions to improve performance of the code

%% Update March 27, 2007 before submission to ICPIG 2007, Prague
%
% # Include N and processes 5, 15, 17 from Benilov and Naidis 2003
% # e+O=e+e+O^+ process 4 from Belilov and Naidis, 2003
% # O^-+NO=NO2 +e detachment process 24 from Benilov and Naidis, 2003
% # Cascade to N2(A) from N2(B) state with account for cascade from N2(C) state to N2(B) state

%% Update November 5, 2006 at EM2C, Paris, France
%
% # Included N2(a') state related processes
% # Included associative ionization processes
% # Included self-quenching of N2(A)
% # Included NO production as a result of N2(A)+O collisions

%% Update October 16, 2005 before URSI GA
%
% # The dissociation rate of O2 was increased by 2.5 in order to match the Naidis 1999 O concentrations

%% Update February 25, 2005
%
% # In this version we expand number of chemical species to model detachment more accurately

%% Update September 22, 2004
%
% # We incude vibrational excitation of N2 molecules and VT relaxation
% processes

%% Update April 20, 2009
%
% # Introduced three missing terms in chemistry model discovered by
%   Jeremy Riousset on April 14, 2009 (kqN2m in eq 7; k15 in eq 11; kdONO
%   in eq 8)
% # Introduced more accurate fast heating and VT heating fractions based
%   on BOLSIG+ calculations by Jeremy Riousset (AirPartition.m) and Ningyu Liu (etaTV.m)
%   (fast heating nut now includes elastic, rotational and 30% of electronic
%   excitation see Naidis 1999 and Popov 2001)
% # Explicitly introduced VT relaxation of N2(v) due to quenching by
%   O using rates of Taylor, Can. J. Chem, 52, 1436, 1974 (see Popov 2001 and Aleksandrov et al
%   J. Phys D Appl Phys. 30 1616, 1997; Aleksandrov et al., Plasma Physics Reports,
%   27, 875, 2001)

%% Update April 25, 2009
%
% # Introduced O2+, O4+ and O2+N2 ions using rates of Kossyi et al.,
%   PSST 1, 207, 1992. The k167 (O2+ + O2 + O2 -> O4+ + O2) process produces O4+
%   and k168 (O2+ + N2 + N2 -> O2+N2 + N2) process produces O2+N2 very fast.
%   O2+N2 is quickly converted to O4+ via k232 (O2+N2 + O2 -> O4+  +
%   N2). O4+ dominates and has a very high recombination rate with electrons
%   (see Alexandrov and Bazelyan PSST 8 285 1999 fig 1 and Kossiy et al)
%   leading to a very significanty delay in streamer to spark
%   transition. Results similar to those of Naidis, J Phys. D Appl. Phys., 32,
%   2649, 1999 can only be reproduced by this version of the model if
%   the above coefficients are set to zero (i.e. assuming O2+ as the only
%   ion, or allowing it to convert to O2+N2 without subsequent conversion into O4+).

%% Update April 30, 2009
%
% # Introduced stopping conditions for ode45 solver when Tg>5000K
%   (sparking condition).
% # Included automatic calculations of t_br for a vector of applied
%   voltages and/or altitudes.

%% Update May 4, 2009
%
% # Corrected value of kqN2A=(1.6+0.77)*1e-10 -> (1.6+0.77)*1e-10 %cm3/s.
%   Summing processes (3) and (4) % of Popov (2001) (repeated hereafter for the sake of clarity),
%   the rate constant should be 2.371e-10 cm3/s in close agreement with Kossyi's second expression
%   for k104 with T = 300 K. Also we note that N2(B) and N2(C) are only excited by electron
%   impact in the current version of the model and converted to N2(A) assuming steady state de-excitation using
%   Einstein's coefficients.
%
% * N2(A) + N2(A) -> N2(B) + N2 (Popov:2001 (3))
% * N2(A) + N2(A) -> N2(C) + N2 (Popov:2001 (4))

%% Update May 5, 2009
%
% # Introduced e-ion recombination for O2+ through the following two
% reactions
%
% * e + O2+ + O2 -> O2 + O2  (Kossyi (44, M=O2))
% * e + O2+ + N2 -> O2 + N2  (Kossyi (44, M=N2))
% * e + O   + O2 -> O  + O2- (Kossyi (48))

%% Update May 6, 2009
%
% # Noticed that a 25% reduction of no in the model provides almost perfect
% agreement with experimental data

%% Update June 1, 2009
%
% # Corrected equations for Kossyi processes I, II, V, VI

%% Update June 2, 2009
%
% # Tried to repeat Dr. Bourdon's results by removing certain reactions

%% Update June 4, 2009
%
% # Change Tn to ni(4) in equation of some constants.

%% Update June 8, 2009 (1)
%
% # added reaction kvi from Kossyi: O2+ + O- + M => O3 + M
% # removed reactions kv: O- + O4+ + M => O + O4 + M (which was wrongly
% interpreted as: O- + O4+ + M => O + O2 + O2 + M)
% # removed reactions kv: O2- + O4+ + M => O2 + O4 + M (which was wrongly
% interpreted as: O2- + O4+ + M => O2 + O2 + O2 + M)
% # removed reactions kv: O2- + O4+ + M => O2 + O4 + M (which was wrongly
% interpreted as: O2- + O4+ + M => O2 + O2 + O2)
% # modified reactions kv: O2+ + O- + M => O + O2 + M
% # included dissociative recombination (ii): O- + O2+ => O + O + O
% # included dissociative recombination (ii): O- + O4+ => O + O2 + O2
% # included dissociative recombination (ii): O- + O2+N2 => O + O2 + N2
% # included dissociative recombination (ii): O2- + O2+ => O2 + O + O
% # included dissociative recombination (ii): O3- + O2+ => O3 + O + O
% # included associative  recombination (i) : O- + O2+ => O + O2
% # added missing 2 factor fo k26ass in eq. (10) and kqN2A in eq. (6)
% # removed formation of O+ and NO+ in eq. (2). Kept formation of N4+ as it
% is readily converted in O2+
% # revised and corrected all rate coefficient rates with addition of a
% reference and number originally associated with the process in the
% reference.

%% Update June 8, 2009 (2)
%
% # Changed calculation of the electron temperature to use Kossyi:1992 and
% air1.m.

%% Update June 10, 2009
%
% # Noticed that k27 and k235 refers to Benilov (21) and Kossyi (235),
% which are the same equation. Removed k235 after comparison of the
% constant rates at different EN_Td and T_K
% # Used Benilov-21 reaction rate instead of Kossyi-57 after sensitivity
% tests. Benilov, Mnatsakanyan and Aleksandrov 1997 all agree. Comtois and
% Kossyi are farther from these values
% # Modified calculation of Fion

%% Update June 12, 2009
%
% # Modified rate for: O2(a) + O2 => O2 + O2 from Lowke:1992 to
% Kossyi:1992-123
% # Included Kossyi-102 in derivation of n1P and n2P

%% Update Oct 12, 2009
%
% # Introduction of parameters Refine and InitialStep in ODE solver
% # 'InitialStep',1e-60 diminishes the initial d.t to a given value (e.g., 1e-60)
% # 'Refine',1e3 forces to multiply the number of steps in the algorithm to
% increase precision
% # Replacememt of the solver ode45 by ode113 to increase precision
% # The above modifications permit to obtain the value t_br=.836mus at
% E=1.9e6V/m when etaT = 1, etaV=0 and no chemistry is involved
% # Correction of the calculation of T_br. The equation is ran until
% Tg>=Tbr. Then evaluation goes through every time step and tbr is found
% when Tg>=Tbr. Tbr is now a parameter passed in integrate6.m.

%% Update Oct 13, 2009 (1)
%
% # Test of eta module
% # Correction of role of QVT. The transfer of vibrational energy of
% Nitrogen to the Translational energy should only apply to N2 and not as
% previously described in v2.9.1. Consequently:
% # QVT is multiplied by nn2/n in equation (4)
% # nn2 replace n in denominator of equation (5)

%% Update Oct 13, 2009 (2)
%
% # Initial test on chemistry. Following discussion with Anne about
% reaction M + e -> O2+ + e + e. She suggested that M can only be equal to
% O2, while we assume that N2+ is also produced and readily converted
% into O2+. Hence her reaction rate is dO2+/dt = kiO2*O2 e, while we have
% dO2+/dt = kiO2 O2 e + kiN2 N2 e. For the sake of comparison we use the
% ionization rate nui for O2 given in Aleksandrov:1995.

%% Update Oct 30, 2009
%
% # Simplification of the calling of each variables
% # kxxx-> k.rxxx

%% Update Nov 23, 2009
%
% # Enclosing the scriptin the function Heating0D.m
% # Uniformisation of the notation between Heating0D.m and Heating1D.m
% # Parameterization for n=cst or p=cst

%% Update Dec 14, 2009
%
% # Fixing of video recording feature using flag.video = 'on'

%% Update Dec 28, 2009
%
% # Introduction of complex ions flags and 3 body chemistry flags.

%% Update Jan 5, 2010
%
% # Reintroduction of the calculation of the rates for all the chemical processes

function [t_br]=Heating1D(U,h,flg)
% Input U (kV), h (km) / Output t_br (s)

global d e eps kB m mu q0 rs x

%% Parameters
N.r      = 501; % Nb of grid points
L        = 2e-3/satm(h); % _m, size of the simulation domain
d.r      = L/(N.r-1); % _m, space step
r        = (0:N.r-1)'*d.r; %_m, domain (assume axisymmetric)
tf       = 1e3; % _s, duration of the simulation
rs       = .2e-3/satm(h); % _m, streamer head radius
T.g0     = 273.15+20; % _K, initial gas temperature (assume isothermal atmosphere)
T.br     = 5000; % _K, breakdown temperature
T.v0     = 273; % _K, initial vibrational temperature of nitrogen
alpha    = .5; % _ , parameter defining the time step
E        = (U-.2)*1e5*satm(h); % _V/_m, applied electric field
x.N2     = .79; % _, fraction of N2(X) in an air mixture
x.O2     = .21; % _, fraction of O2(X) in an air mixture
eps      = 1e-7; % _, relative tolerance in convergence criterion 

n.e0     = 2e20*satm(h)*satm(h); % _m-3, initial streamer electron density
n.O2p0   = n.e0; % _m-3, initial density of O2+
n.O4p0   = 0; % _m-3, initial density of O4+
n.O2pN20 = 0; % _m-3, initial density of O2+N2
n.Om0    = 0; % _m-3, initial density of O-
n.O2m0   = 0; % _m-3, initial density of O2-
n.O2a0   = 0; % _m-3, initial density of O2(a)
n.N2a0   = 0; % _m-3, initial density of N2(a')
n.N2A0   = 0; % _m-3, initial density of N2(A)
n.O0     = 0; % _m-3, 2e10 at 70 km initial density of O
n.N0     = 0; % _m-3, initial density of N
n.NO0    = 0; % _m-3, initial density of NO: 2e8 at 70 km; 
n.N2p0   = 0; % _m-3, initial density of N2+
n.N4p0   = 0; % _m-3, initial density of N4+
n.O30    = 0; % _m-3, initial density of O3

n.out    = []; % [_t, (m^-3)^12], output matrix with species densities on the axis
global FF
FF.iz    = []; % _s^-1, ionization rate
FF.st    = []; % _s^-1, stepwise ionisation rate
FF.a2    = []; % _s^-1, 2-body attachment rate
FF.a3    = []; % _s^-1, 3-body attachment rate
FF.ei    = []; % _s^-1, electron-ion recombination rate
FF.de    = []; % _s^-1, detachment rate
FF.ex    = []; % _s^-1, electron impact excitation rate
FF.di    = []; % _s^-1, electron impact dissociation rate
FF.gd    = []; % _s^-1, ground states chemistry rate
FF.rd    = []; % _s^-1, radicals chemistry rate
FF.ii    = []; % _s^-1, ion-ion recombination rate
FF.ip    = []; % _s^-1, positive ion chemistry rate
FF.im    = []; % _s^-1, negative ion chemistry rate

%% Physical constants
gamma    = 1.4; % _, specific heat ratio of air (assume perfect gas)
kB       = 1.38e-23; % _J/_K, Boltzmann constant
e        = 1.6e-19; % _C, electron charge
m.g      = 4.82e-26; % _kg, average mass of air molecule
p.gnd    = 1e5; % _Pa, unperturbed pressure at ground level
d.Ev     = 0.29*e; % _J, vibrational quanta
Ev.g0    = 1.291e-6; % _eV, starting mean vibrational energy (a small number), Tv = 300 K

n.gnd    = 2.5e25;%p.gnd/kB/T.g0; % _m^-3, unperturbed density at ground level

%% Calculated parameters
n.g0     = n.gnd*satm(h); % _m^-3, unperturbed density at altitude h 
n.N20    = x.N2*n.g0; % _m-3, molecular nitrogen density
n.O20    = x.O2*n.g0; % _m-3, molecular oxygen density
rho.g0   = n.g0*m.g; %1.2041*satm(h); %_kg/m^3, unpertubed mass density
v.g0     = 0; % _m/s, unperturbed molecules velocity
p.g0     = n.g0*kB.*T.g0; %101325*satm(h); %_Pa, atmospheric pressure
Ev.g0    = d.Ev/(exp(d.Ev/(kB*T.v0))-1); % _J, starting mean vibrational energy (a small number) [Naidis:2007]

cs.o     = sqrt(gamma*p.g0/rho.g0); %_m/_s, speed of sound
d.Ev_eV  = d.Ev/e; % _eV, vibrational quanta

%% Variables
t        = 0; % _s, current time
rho.g    = ones(N.r,1)*rho.g0; %_kg/m3, current gas mass density
v.g      = ones(N.r,1)*v.g0; % _v, current neutral particules velocity
p.g      = ones(N.r,1)*p.g0; % _Pa, current total pressure (obtained from n.g, T.g)
Ev.g     = ones(N.r,1)*Ev.g0; % _J, Energy transfered to the vibrational levels of nitrogen
T.g      = ones(N.r,1)*T.g0; % _K, current gas temperature

nn = 1;
set(gcf,'Units','normalized','OuterPosition',[1/3 1/3 2/3 2/3],'Color',[1 1 1])
while(t<=tf && max(T.g)<= T.br)
    cs.max   = max((gamma*p.g./rho.g).^.5);
    d.t      = alpha*min([d.r/max(cs.o,cs.max),min(tauVT(rho.g,p.g))]);%,1e-9]);
    
    % On-axis values for use in accoustic model, index "0" indicates this
    n.g0     = rho.g(1)/m.g; %_m^-3, (needs converted to _cm^-3 before used in the chemistry model)
    T.g0     = p.g(1)/(n.g0*kB); % _K 
    Ev.g0    = Ev.g(1); % _J
    % mu.e     = air(Ef,n.g0*1e-6,11); % _V_m2/_s, electron mobility (Zhao:1992)
    mu.e     = morrowair(E,n.g0*1e-6,4); %_V_m2/_s, electron mobility (Zhao:1992)
    mu.p     = 2.5e-4/(n.g0/n.gnd); % _V_m2/_s, positive ion mobility (Zhao:1992)
    mu.n     = 2.2e-4/(n.g0/n.gnd); % _V_m2/_s, negative ion mobility (Zhao:1992)

    % Run chemistry during each interval [t,t+d.t] assuming n.g, T.g, and
    % Ev.g are constant on this interval of time
    u0       = [n.e0;n.O2p0;n.O4p0;n.O2pN20;n.Om0;n.O2m0;n.O2a0;n.N2a0;n.N2A0;n.O0;n.N0;n.NO0;n.N2p0;n.N4p0;n.O30]; % _m3

    [tau,u]  = chemistry1(t, t+d.t, n.gnd*1e-6, E, u0*1e-6, n.g0*1e-6, T.g0, Ev.g0, flg); % _s, _cm3, output all densities in cm3!!!
    n.e0     = u(end,1)*1e6; % _m^-3
    n.O2p0   = u(end,2)*1e6; % _m^-3
    n.O4p0   = u(end,3)*1e6; % _m^-3
    n.O2pN20 = u(end,4)*1e6; % _m^-3
    n.Om0    = u(end,5)*1e6; % _m^-3
    n.O2m0   = u(end,6)*1e6; % _m^-3
    n.O2a0   = u(end,7)*1e6; % _m^-3
    n.N2A0   = u(end,8)*1e6; % _m^-3
    n.N2a0   = u(end,9)*1e6; % _m^-3
    n.O0     = u(end,10)*1e6; % _m^-3
    n.N0     = u(end,11)*1e6; % _m^-3
    n.NO0    = u(end,12)*1e6; % _m^-3
    n.N2p0   = u(end,13)*1e6; % _m^-3
    n.N4p0   = u(end,14)*1e6; % _m^-3
    n.O30    = u(end,10)*1e6; % _m^-3
    n.O3m0   = max(n.O2p0+n.O4p0+n.O2pN20+n.N2p0+n.N4p0-n.e0-n.Om0-n.O2m0, 0); % _m^-3
    
    % Calculate values (average of the values at t and t+d.t)
    n.ei     = (u(1,1)+u(end,1))/2*1e6; % _m^-3
    n.O2pi   = (u(1,2)+u(end,2))/2*1e6; % _m^-3
    n.O4pi   = (u(1,3)+u(end,3))/2*1e6; % _m^-3
    n.O2pN2i = (u(1,4)+u(end,4))/2*1e6; % _m^-3
    n.Omi    = (u(1,5)+u(end,5))/2*1e6; % _m^-3
    n.O2mi   = (u(1,6)+u(end,6))/2*1e6; % _m^-3
    n.O2ai   = (u(1,7)+u(end,7))/2*1e6; % _m^-3
    n.N2ai   = (u(1,8)+u(end,8))/2*1e6; % _m^-3
    n.N2Ai   = (u(1,9)+u(end,9))/2*1e6; % _m^-3
    n.Oi     = (u(1,10)+u(end,10))/2*1e6; % _m^-3
    n.Ni     = (u(1,11)+u(end,11))/2*1e6; % _m^-3
    n.NOi    = (u(1,12)+u(end,12))/2*1e6; % _m^-3
    n.N2pi   = (u(1,13)+u(end,13))/2*1e6; % _m^-3
    n.N4pi   = (u(1,14)+u(end,14))/2*1e6; % _m^-3
    n.O3i    = (u(1,15)+u(end,15))/2*1e6; % _m^-3
    n.O3mi   = max(n.O2pi+n.O4pi+n.O2pN2i+n.N2pi+n.N4pi-n.ei-n.Omi-n.O2mi, 0); %this avoids negative nO3m
    
    sigma    = e*(mu.p*(n.O2pi+n.O4pi+n.O2pN2i+n.N2pi+n.N4pi)+mu.e*n.ei+mu.n*(n.O2pi+n.O4pi+n.O2pN2i-n.ei)); %S/m, on axis
    I        = sigma*E*pi*rs^2;
    q0       = sigma*E^2; % _J/m^3/s, Joule heating on axis, q(r) = q.o*exp(-r^2/rs^2);
    
    % 1-D values for LW update
    eta.T    = AirPartition(E,n.g0*1e-6,15); %_, fraction of Joule energy transfered to fast heating (Based on Bolsig+)
    eta.V    = AirPartition(E,n.g0*1e-6,7); %_, fraction of Joule energy transfered to vibrational level of N2 (Based on Bolsig+)
    
    dTdtFH   = eta.T*(gamma-1)*q0/(n.g0*kB);
    dTdtVT   =       (gamma-1)*qVT(rho.g0,p.g0,Ev.g0)/(n.g0*kB);
    
    n.out    = [n.out; t+d.t/2 n.ei n.O2pi n.O4pi n.O2pN2i n.Omi n.O2mi n.O3mi n.O2ai n.N2ai n.N2Ai n.Oi n.Ni n.NOi n.N2pi n.N4pi n.O3i n.g0 T.g0 d.Ev/(kB*log(1+d.Ev/Ev.g0)) I dTdtFH dTdtVT];
        
    [t,rho,v,p,Ev]=LW(t,rho,v,p,Ev);
    t  = t + d.t;
    
    T.g = p.g*m.g./rho.g/kB;

    if(strcmp(flg.plot,'on'))
        clf;
        fprintf('t = %2.3es\n',t);
        % subplot(331);
        % plot(r,rho.g,'k')
        % axis([r(1) r(N.r) .75 1.5])
        % box on
        % ylabel('\rho (kg/m^3)');
        % title(['t = ',num2str(t,'%2.2e'),' s'])

        subplot(331);
        plot(r,rho.g/rho.g0,'k')
        axis([r(1) r(N.r) .75 1.1])
        box on
        ylabel('n_g/n_{g_0}');
        title(['t = ',num2str(t,'%2.2e'),' s'])
        set(gca,'XMinorTick','on','YMinorTick','on','TickDir','out')
        
        subplot(334);
        plot(r,v.g,'k-')
        axis([r(1) r(N.r) 0 125])
        ylabel('v (m/s)');
        set(gca,'XMinorTick','on','YMinorTick','on','TickDir','out')
        
        subplot(337);
        semilogy(r,p.g,'k-')
        axis([r(1) r(N.r) 1 20*p.gnd])
        box on
        xlabel('r (m)')
        ylabel('p (Pa)');
        set(gca,'XMinorTick','on','YMinorTick','on','TickDir','out')
        
        subplot(3,3,[2,5,8]);
        semilogy(r,T.g,'k-')
        axis([r(1) r(N.r) 0 T.br])
        box on
        xlabel('r (m)')
        ylabel('T (K)');
        set(gca,'XMinorTick','on','YMinorTick','on','TickDir','out')
        
        subplot(3,3,[3,6,9]);
        semilogy(r,Ev.g/e,'k-')
        axis([r(1) r(N.r) 1e-6 1e1])
        box on
        xlabel('r (m)')
        ylabel('E_v (eV)');
        set(gca,'XMinorTick','on','YMinorTick','on','TickDir','out')
        
        if(strcmp(flg.movie,'on'))
            M(nn) = getframe(gcf);
        else
            getframe(gcf);
        end
        nn    = nn+1;
    end
end
if(strcmp(flg.movie,'on'))
    movie2avi(M,[num2str(h),'km.avi'],'quality',100);
end
   
t_br = t(end);
if(strcmp(flg.plot,'on'))
    t       = n.out(:,1);
    n.e     = n.out(:,2);

    n.O2p   = n.out(:,3);
    n.O4p   = n.out(:,4);
    n.O2pN2 = n.out(:,5);
    n.N2p   = n.out(:,15);
    n.N4p   = n.out(:,16);
    
    n.Om    = n.out(:,6);
    n.O2m   = n.out(:,7);
    n.O3m   = n.out(:,8);

    n.O2a   = n.out(:,9);
    n.N2A   = n.out(:,10);
    n.N2a   = n.out(:,11);

    n.O     = n.out(:,12);
    n.N     = n.out(:,13);
    n.NO    = n.out(:,14);
    n.O3    = n.out(:,17);
    
    n.gout  = n.out(:,18);
    T.gout  = n.out(:,19);
    T.vout  = n.out(:,20);
    I       = n.out(:,21);
    
    dT_dtFH = n.out(:,22);
    dT_dtVT = n.out(:,23);
    
    % Plots
    figure(2)
    clf(2)
    set(gcf,'Color',[1 1 1])
    set(gcf,'Units','Normalized','OuterPosition',[0 0 1 1])
    subplot(3,3,[1,4,7])
    loglog(t,n.e,'-',t,n.O2p,'+-', t,n.O4p,'x-',t,n.O2pN2,'*-',t,n.N2p,'p-',t,n.N4p,'h-',t,n.Om,'--', t,n.O2m,'-.', t,n.O3m,':')
    xlim([t_br/1e3 t_br] )
    %ylim([1e9 1e15]*satm(h))
    xlabel('Time (s)', 'FontSize', 12)
    ylabel('Number density (cm^{-3})', 'FontSize', 12)
    legend('e', 'O_2^+', 'O_4^+', 'O_2^+N_2', 'N_2^+', 'N_4^+', 'O^-', 'O_2^-', 'O_3^-', 0);
    set(gca, 'FontSize', 12,'XminorTick','on','YminorTick','on')
    legend('boxoff');
    title('Electrons and Ions')
    % grid

    subplot(3,3,[2,5,8])
    loglog(t,FF.iz,'-', t,FF.st,'--', t,FF.a2+FF.a3,'-.',t,FF.ei,':', t,FF.de,'+-', t,FF.ex,'x-', t,FF.di,'p-', t, FF.gd, 'h-', t,FF.ii,'*-', t,FF.ip,'d-', t,FF.im,'s-', t,FF.rd,'o-');
    xlim([t_br/1e3 t_br] )
    %ylim([1e4 1e10]*satm(h))
    xlabel('Time (s)', 'FontSize', 12)
    ylabel('Rates (1/s)', 'FontSize', 12)
    legend('direct iz', 'stepwise iz', 'att', 'e-i rec', 'detachment', 'radicals excitation', 'dissociation', 'gnd states chem', 'i-i rec', 'M^+ chem', 'M^- chem', 'rad chem', 0);
    legend('boxoff');
    set(gca,'FontSize', 12,'XminorTick','on','YminorTick','on')
    % grid
    
    subplot(3,3,3)
    loglog(t,n.O2a,'-', t,n.N2a,'-.', t, n.O,'-',t,n.N2A,'--', t, n.NO, '-',t,n.N,'--',t,n.O3,'k:')
    xlim([t_br/1e3 t_br] )
    %ylim([1e11 1e17]*satm(h))
    xlabel('Time (s)', 'FontSize', 12)
    ylabel('Number density (cm^{-3})', 'FontSize', 12)
    legend('O_2(a)', 'N_2(a`)', 'O', 'N_2(A)', 'NO', 'N', 'O_3', 0);
    legend('boxoff');
    set(gca, 'FontSize', 12,'XminorTick','on','YminorTick','on')
    title('Neutral and Metastable Species')
    grid

    subplot(3,3,6)
    semilogx(t,T.gout,'-', t, T.vout, '--')
    xlim([t_br/1e3 t_br] )
    ylim([T.g0 T.br])
    xlabel('Time (s)', 'FontSize', 12)
    ylabel('Temperature (K)', 'FontSize', 12)
    legend('T_g', 'T_v', 0);
    legend('boxoff');
    set(gca,'FontSize', 12,'XminorTick','on','YminorTick','on')
    grid
    
    subplot(3,3,9)
    loglog(t,I)
    xlim([t_br/1e3 t_br] )
    xlabel('Time (s)', 'FontSize', 12)
    ylabel('I (A)', 'FontSize', 12)
    legend('boxoff');
    set(gca, 'FontSize', 12,'XminorTick','on','YminorTick','on')
    grid

    figure(3)
    clf(3)
    set(gcf,'Color',[1 1 1])
    set(gcf,'Units','Normalized','OuterPosition',[0 0 1 1])
    subplot(3,1,1)
    semilogx(t,n.gout*kB.*T.gout*1e-5,'k')
    xlim([t_br/1e3 t_br] )
    xlabel('Time (s)', 'FontSize', 12)
    ylabel('Pressure (bar)', 'FontSize', 12)
    legend('boxoff');
    set(gca,'FontSize', 12,'XminorTick','on','YminorTick','on')
    grid

    subplot(3,1,2)
    semilogx(t,n.gout/n.gnd,'-')
    xlim([t_br/1e3 t_br] )
    xlabel('Time (s)', 'FontSize', 12)
    ylabel('n_g/n_{gnd}', 'FontSize', 12)
    legend('boxoff');
    set(gca,'FontSize', 12,'XminorTick','on','YminorTick','on')
    grid
    
    subplot(3,1,3)
    loglog(t,dT_dtFH,'r',t,dT_dtVT,'b--')
    xlim([t_br/1e3 t_br] )
    xlabel('Time (s)', 'FontSize', 12)
    ylabel('dT/dt (K/s)', 'FontSize', 12)
    legend( 'Fast heating rate', 'VT heating rate');
    legend('boxoff');
    set(gca,'FontSize', 12,'XminorTick','on','YminorTick','on')
    grid

end

fprintf('\n\tt_br = %2.2e s.\n\n',t_br)

    function [t,rho,v,p,Ev]=LW(t,rho,v,p,Ev)
        % * Lax Wendroff (2 steps)
        % * Non-linear, cylindrical geometry
        % * Using a RHS side (ie, decomposing the divergence term into a Cartesian-
        % like term and a RHS).

        old.rho       = rho.g; %_kg/_m^3
        old.v         = v.g; % _m/_s
        old.p         = p.g; % _Pa
        old.Ev        = Ev.g; % _J

        hf1.rho       = zeros(N.r,1); % _kg/_m^3
        hf1.v         = zeros(N.r,1); % _m/_s
        hf1.p         = zeros(N.r,1); % _Pa
        hf1.Ev        = zeros(N.r,1); % _J
        hf2.rho       = zeros(N.r,1); % _m^-3
        hf2.v         = zeros(N.r,1); % _m/_s
        hf2.p         = zeros(N.r,1); % _Pa
        hf2.Ev        = zeros(N.r,1); % _J
        F.in.rho      = zeros(N.r,1); % _kg/_m^3
        F.in.v        = zeros(N.r,1); % _m/_s
        F.in.p        = zeros(N.r,1); % _Pa
        F.in.Ev       = zeros(N.r,1); % _J
        F.out.rho     = zeros(N.r,1); % _kg/_m^3
        F.out.v       = zeros(N.r,1); % _m/_s
        F.out.p       = zeros(N.r,1); % _Pa
        F.out.Ev      = zeros(N.r,1); % _J
        RHS.rho       = zeros(N.r,1); % _kg/_m3
        RHS.v         = zeros(N.r,1); % _m/_s
        RHS.p         = zeros(N.r,1); % _Pa
        RHS.Ev        = zeros(N.r,1); % _J

        ii            = 2:N.r-1;
    
        %% LW non-linear
        rho.g(1) = old.rho(2);
        v.g(1)   = old.v(2);
        p.g(1)   = old.p(2);
        Ev.g(1)  = old.Ev(2);

        F.in.rho(ii)  = old.rho(ii-1).*old.v(ii-1);
        F.out.rho(ii) = old.rho(ii)  .*old.v(ii);
        RHS.rho(ii)   = (old.rho(ii).*old.v(ii)+old.rho(ii-1).*old.v(ii-1))/2./(r(ii)-d.r/2);
        hf1.rho(ii)   = (old.rho(ii) + old.rho(ii-1))/2 ...
            - d.t/d.r/2*(F.out.rho(ii) - F.in.rho(ii) ) ...
            - d.t/2*RHS.rho(ii) ;

        F.in.rho(ii)  = old.rho(ii)  .*old.v(ii);
        F.out.rho(ii) = old.rho(ii+1).*old.v(ii+1);
        RHS.rho(ii)   = (old.rho(ii).*old.v(ii)+old.rho(ii+1).*old.v(ii+1))/2./(r(ii)+d.r/2);
        hf2.rho(ii)   = (old.rho(ii) + old.rho(ii+1))/2 ...
            - d.t/d.r/2*(F.out.rho(ii) - F.in.rho(ii)) ...
            - d.t/2*RHS.rho(ii) ;
        
        F.in.v(ii)    = (old.v(ii-1)).^2/2;
        F.out.v(ii)   = (old.v(ii))  .^2/2;
        RHS.v(ii)     = (old.p(ii)-old.p(ii-1))/d.r./((old.rho(ii)+old.rho(ii-1))/2);
        hf1.v(ii)     = (old.v(ii) + old.v(ii-1))/2 ...
            - d.t/d.r/2*(F.out.v(ii) - F.in.v(ii)) ...
            - d.t/2*RHS.v(ii) ;

        F.in.v(ii)    = (old.v(ii))  .^2/2;
        F.out.v(ii)   = (old.v(ii+1)).^2/2;
        RHS.v(ii)     = (old.p(ii+1)-old.p(ii))/d.r./((old.rho(ii+1)+old.rho(ii))/2);
        hf2.v(ii)     = (old.v(ii+1) + old.v(ii))/2 ...
            - d.t/d.r/2*(F.out.v(ii) - F.in.v(ii)) ...
            - d.t/2*RHS.v(ii) ;

        F.in.p(ii)    = gamma*old.p(ii-1).*old.v(ii-1);
        F.out.p(ii)   = gamma*old.p(ii)  .*old.v(ii);
        RHS.p(ii)     = gamma*(old.p(ii-1).*old.v(ii-1)+old.p(ii).*old.v(ii))/2./(r(ii)-d.r/2) ...
            - (gamma-1)*((old.v(ii-1)+old.v(ii))/2.*(old.p(ii)-old.p(ii-1))/d.r ...
            + eta.T.*(q(r(ii-1))+q(r(ii)))/2 ...
            + (qVT(old.rho(ii-1),old.p(ii-1),old.Ev(ii-1))+qVT(old.rho(ii),old.p(ii),old.Ev(ii)))/2);
        hf1.p(ii)     = (old.p(ii) + old.p(ii-1))/2 ...
            - d.t/d.r/2*(F.out.p(ii) - F.in.p(ii) ) ...
            - d.t/2*RHS.p(ii) ;

        F.in.p(ii)    = gamma*old.p(ii)  .*old.v(ii);
        F.out.p(ii)   = gamma*old.p(ii+1).*old.v(ii+1);
        RHS.p(ii)     = gamma*(old.p(ii).*old.v(ii)+old.p(ii+1).*old.v(ii+1))/2./(r(ii)+d.r/2) ...
            - (gamma-1)*((old.v(ii)+old.v(ii+1))/2.*(old.p(ii+1)-old.p(ii))/d.r ...
            + eta.T.*(q(r(ii+1))+q(r(ii)))/2 ...
            + (qVT(old.rho(ii+1),old.p(ii+1),old.Ev(ii+1))+qVT(old.rho(ii),old.p(ii),old.Ev(ii)))/2);
        hf2.p(ii)     = (old.p(ii+1) + old.p(ii))/2 ...
            - d.t/d.r/2*(F.out.p(ii) - F.in.p(ii) ) ...
            - d.t/2*RHS.p(ii) ;

        F.in.Ev(ii)   = old.Ev(ii-1).*old.v(ii-1);
        F.out.Ev(ii)  = old.Ev(ii)  .*old.v(ii);
        RHS.Ev(ii)    = (old.Ev(ii-1).*old.v(ii-1)+old.Ev(ii).*old.v(ii))/2./(r(ii)-d.r/2) ...
            - eta.V.*(q(r(ii-1))+q(r(ii)))/2./(x.N2*(old.rho(ii)+old.rho(ii-1))/2/m.g) ...
            + (qVT(old.rho(ii-1),old.p(ii-1),old.Ev(ii-1))+qVT(old.rho(ii),old.p(ii),old.Ev(ii)))/2./(x.N2*(old.rho(ii)+old.rho(ii-1))/2/m.g);
        hf1.Ev(ii)    = (old.Ev(ii) + old.Ev(ii-1))/2 ...
            - d.t/d.r/2*(F.out.Ev(ii) - F.in.Ev(ii) ) ...
            - d.t/2*RHS.Ev(ii) ;

        F.in.Ev(ii)   = old.Ev(ii)  .*old.v(ii);
        F.out.Ev(ii)  = old.Ev(ii+1).*old.v(ii+1);
        RHS.Ev(ii)    = (old.Ev(ii).*old.v(ii)+old.Ev(ii+1).*old.v(ii+1))/2./(r(ii)+d.r/2) ...
            - eta.V.*(q(r(ii))+q(r(ii+1)))/2./(x.N2*(old.rho(ii)+old.rho(ii+1))/2/m.g) ...
            + (qVT(old.rho(ii),old.p(ii),old.Ev(ii))+qVT(old.rho(ii+1),old.p(ii+1),old.Ev(ii+1)))/2./(x.N2*(old.rho(ii)+old.rho(ii+1))/2/m.g);
        hf2.Ev(ii)    = (old.Ev(ii) + old.Ev(ii+1))/2 ...
            - d.t/d.r/2*(F.out.Ev(ii) - F.in.Ev(ii) ) ...
            - d.t/2*RHS.Ev(ii) ;

        F.in.rho(ii)  = hf1.rho(ii).*hf1.v(ii);
        F.out.rho(ii) = hf2.rho(ii).*hf2.v(ii);
        RHS.rho(ii)   = (hf1.rho(ii).*hf1.v(ii)+hf2.rho(ii).*hf2.v(ii))/2./r(ii);
        rho.g(ii)       = (old.rho(ii) ...
            - d.t/d.r*(F.out.rho(ii) - F.in.rho(ii) ) ...
            - d.t*RHS.rho(ii));
        
        F.in.v(ii)    = (hf1.v(ii)).^2/2;
        F.out.v(ii)   = (hf2.v(ii)).^2/2;
        RHS.v(ii)     = (hf2.p(ii)-hf1.p(ii))/d.r./((hf1.rho(ii)+hf2.rho(ii))/2);
        v.g(ii)       = old.v(ii) + ...
            - d.t/d.r*(F.out.v(ii) - F.in.v(ii)) ...
            - d.t*RHS.v(ii) ;

        F.in.p(ii)    = gamma*hf1.p(ii).*hf1.v(ii);
        F.out.p(ii)   = gamma*hf2.p(ii).*hf2.v(ii);
        RHS.p(ii)     = gamma*(hf1.p(ii).*hf1.v(ii)+hf2.p(ii).*hf2.v(ii))/2./r(ii) ...
            - (gamma-1)*((hf1.v(ii)+hf2.v(ii))/2.*(hf2.p(ii)-hf1.p(ii))/d.r ...
            + eta.T.*(q(r(ii)-d.r/2)+q(r(ii)+d.r/2))/2 ...
            + (qVT(hf1.rho(ii),hf1.p(ii),hf1.Ev(ii))+qVT(hf2.rho(ii),hf2.p(ii),hf2.Ev(ii)))/2);
        p.g(ii)       = (old.p(ii) ...
            - d.t/d.r*(F.out.p(ii) - F.in.p(ii) ) ...
            - d.t*RHS.p(ii)) ;

        F.in.Ev(ii)   = hf1.Ev(ii).*old.v(ii);
        F.out.Ev(ii)  = hf2.Ev(ii).*old.v(ii);
        RHS.Ev(ii)    = (hf1.Ev(ii).*hf1.v(ii)+hf2.Ev(ii).*hf2.v(ii))/2./r(ii) ...
            - eta.V.*(q(r(ii)-d.r/2)+q(r(ii)+d.r/2))/2./(x.N2*(hf1.rho(ii)+hf2.rho(ii))/2/m.g) ...
            + (qVT(hf1.rho(ii),hf1.p(ii),hf1.Ev(ii))+qVT(hf2.rho(ii),hf2.p(ii),hf2.Ev(ii)))/2./(x.N2*(hf1.rho(ii)+hf2.rho(ii))/2/m.g);
        Ev.g(ii)      = old.Ev(ii) ...
            - d.t/d.r*(F.out.Ev(ii) - F.in.Ev(ii) ) ...
            - d.t*RHS.Ev(ii) ;

        rho.g(N.r)    = old.rho(N.r-1);
        v.g(N.r)      = old.v(N.r-1);
        p.g(N.r)      = old.p(N.r-1);
        Ev.g(N.r)     = old.Ev(N.r-1);
    end
end