function r=alexair(e,h,i)

%	function r=alexair(e,h,i)
%
%	Calculates ionization coefficients of N2 and O2, two body attachment
%       coefficient, and various excitation coefficeints in air
%       at different applied electric fields.
%       
%	r - ionization/attacment/excitation 1/s;
%	e - array electric field in V/m (for valid results 0<e*no/n<2.68e7 V/m, no=2.688e19 cm-3, or E/N<1000 Td);
%	h - altitude in km (h>150 km is intepreted as neutral density n in cm-3);
%	i - integer parameter specifying process to be calculated:
%
%       1 Ionization N2				
%       2 Ionization O2		
%       3 Two-body attachment					
%       4 A^3 Sigma_u^+ N2
%       5 B^3 Pi_g N2 (1st positive N2)
%       6 a^1 Pi_g N2
%       7 C^3 Pi_u N2 (2nd positive N2)
%       8 a^1 Delta_g O2
%       9 b^1 Sigma_g^+ O2
%       10 O2 Dissociation 
%       11 O2 Dissociation with O(^1D) production 
%       12 a'^1 Sigma_u^- N2 
%
%	The function is created by Victor Pasko (CSSL Lab,
%	Penn State University) on July 6, 2005 and is based on 
%       analytical functions specified in 
%       Alexandrov, N. L., A. E. Bazelyan, E. M. Bazelyan, and I. V. Kochetov,
%       Modeling of long streamers in atmospheric-pressure air,
%       Plasma Phys. Rep., 21, 57-75, 1995.
%
%       Updated on November 3, 2006 by adding N2(a'^1 Sigma_u^-) 
%       using approach detailed in Kossyi et al., PSST, 1, 207, 1992  
%
%Define some constants:
N0=2.688e19; %neutral density cm-3 at the ground at temperature 273 deg. K

%Determine neutral density N in cm-3:

if(0<=h && h<150), %means input is altitude in km

%US Standard Atmosphere:
%We replaced original 2.5e19cm-3 at 0 km altitude
%by 2.688e19 cm-3 which is our reference
%number density at temperature 273 K.
%statm(1,:) altitude in km
%statm(2,:) neutral density in cm-3

statm=[
 0e+0	 2.688e+19
 5e+0	 1.53e+19
 1e+1	 8.59e+18
 1.5e+1	 4.05e+18
 2e+1	 1.85e+18
 2.5e+1	 8.33e+17
 3e+1	 3.83e+17
 3.5e+1	 1.76e+17
 4e+1	 8.31e+16
 4.5e+1	 4.088e+16
 5e+1	 2.13e+16
 5.5e+1	 1.181e+16
 6e+1	 6.33e+15
 6.4e+1	 3.93e+15
 6.8e+1	 2.39e+15
 7.2e+1	 1.39e+15
 7.6e+1	 7.72e+14
 8e+1	 4.03e+14
 8.4e+1	 1.99e+14
 8.8e+1	 9.48e+13
 9.2e+1	 4.37e+13
 9.6e+1	 2.07e+13
 1e+2	 1.04e+13
 1.08e+2	 3.18e+12
 1.14e+2	 1.43e+12
 1.2e+2	 6.61e+11
 1.26e+2	 3.4e+11
 1.32e+2	 1.91e+11
 1.4e+2	 9.7e+10
 1.5e+2	 4.92e+10
];
	N=10^(interp1(statm(:,1),log10(statm(:,2)),h));
else%means input is neutral density in cm-3
	N=h;
end %if

% E/N in V*cm2
EN=e/N*1e-2; % 1e-2 since e is in V/m
g=EN*1e16; % V*cm2, gamma from Alexandrov paper

%Oxygen and nitrogen molecular densities in cm-3
N2=0.79*N;
O2=0.21*N;

%Define rates

switch i
    case 1 %ionization N2 [1/s]
        r=(g>=8).*(g<=30).*10.^(-8.09-40.29./g)+...
          (g>30).*10.^(-7.37-61.81./g);
        r=r*N2; %1/s
    case 2 %ionization O2 [1/s]
        r=(g>=6).*(g<=26).*10.^(-8.31-28.57./g)+...
          (g>26).*10.^(-7.54-48.57./g+log10(1+4e-7*g.^3));
        r=r*O2; %1/s               
    case 3 %Two-body attachment
        r=(g>=3).*(g<=9).*10.^(-9.42-12.7./g)+...
          (g>9).*(g<=30).*10.^(-10.21-5.7./g)+...
          (g>30).*10.^(-10.33-2.7e-3.*g);
        r=r*O2; %1/s   
    case 4 %A^3 Sigma_u^+
        r=(g>=2).*(g<=10).*10.^(-8.4-17.11./g)+...
          (g>10).*10.^(-8.67-14.41./g);
        r=r*N2; %1/s   
    case 5 %B^3 Pi_g
        r=(g>=2).*(g<=10).*10.^(-7.91-16.81./g)+...
          (g>10).*10.^(-8.2-13.92./g);
        r=r*N2; %1/s         
    case 6 %a^1 Pi_g
        r=(g>=3).*(g<=10).*10.^(-8.17-18.74./g)+...
          (g>10).*10.^(-8.29-17.53./g);
        r=r*N2; %1/s 
    case 7 %C^3 Pi_u
        r=(g>=4).*(g<=20).*10.^(-7.88-23.32./g)+...
          (g>20).*10.^(-8.08-19.37./g);
        r=r*N2; %1/s 
    case 8 %a^1 Delta_g
        r=(g>=2).*(g<=4.5).*10.^(-10.26-0.79./g)+...
          (g>4.5).*(g<=20).*10.^(-8.74-7.68./g)+...
          (g>20).*10.^(-8.68-7.7./g-3.1e-3.*g);
        r=r*O2; %1/s 
    case 9 %b^1 Sigma_g^+
        r=(g>=2).*(g<=5).*10.^(-12.06+3.97e-2.*g.^2)+...
          (g>5).*(g<=10).*10.^(-9.13-9.7./g)+...
          (g>10).*10.^(-9.4-7.05./g-1.7e-3.*g);
        r=r*O2; %1/s 
    case 10 %O2 Dissociation 
        r=(g>=3).*(g<=10).*10.^(-7.78-14.08./g)+...
          (g>10).*(g<=30).*10.^(-8.31-8.78./g)+...
          (g>30).*2.5e-9;
        r=r*O2; %1/s 
    case 11 %O2 Dissociation with O(^1D) production 
        r=(g>=3).*(g<=10).*10.^(-7.43-17.06./g)+...
          (g>10).*10.^(-7.6-15.43./g);
          r=r*O2; %1/s  
    case 12 %a'^1 Sigma_u^-
        r=10.^(-8.8-16.7./g)+...
          10.^(-8.5-17.4./g)+...
          10.^(-8.7-17.5./g);
          r=r*N2; %1/s   

    otherwise
        disp('Wrong input i')
end

if(max(g) > 100)
    fprintf('Warning: ELECTRIC FIELD IS TOO LARGE (E/N = %1.2e Td > 1000 Td), RESULT MAY NOT BE VALID.\n',e/N/100/1e-17)
end %if














