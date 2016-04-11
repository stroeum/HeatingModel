function r=morrowair(e,h,i)

%	function r=morrowair(e,h,i)
%
%	Calculates ionization coefficient, two and three body attachment
%       coefficients, electrtron mobility and electron diffusion in air
%       at different applied electric fields.
%       
%	r - ionization/attacment 1/s, or mobility m^2/V/s, or diffusion m^2/s;
%	e - array electric field in V/m (for valid results 0<e*no/n<6.3e7 V/m, no=2.688e19 cm-3);
%	h - altitude in km (h>150 km is intepreted as neutral density n in cm-3);
%	i - integer parameter specifying process to be calculated:
%
%		1 Ionization				
%		2 Two body attachment		
%		3 Three body attachment					
%		4 Mobility
%		5 Diffusion
%       6 Babaeva 2 body attachment
%       7 Babaeva 3 body attachment
%
%	The function is created by Victor Pasko (CSSL Lab,
%	Penn State University) on August 22, 2003 and is based on 
%       analytical functions specified in 
%       Morrow, R., and J. J. Lowke, Streamer propagation in air, 
%       J. Phys. D: Appl. Phys., 30, 614-627, 1997.


%Define some constants:
N0=2.688e19; %neutral density cm-3 at the ground at temperature 273 deg. K

%Determine neutral density N in cm-3:

if(0<=h & h<150), %means input is altitude in km

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

% Absolute value of electron drift velocity in cm/s
We=(EN>2e-15).*(7.4e21*EN+7.1e6)+...
   (EN>=1e-16).*(EN<=2e-15).*(1.03e22*EN+1.3e6)+...
   (EN>=2.6e-17).*(EN<=1e-16).*(7.2973e21*EN+1.63e6)+...
   (EN<=2.6e-17).*(6.87e22*EN+3.38e4);

%Define rates

switch i
    case 1 %ionization 1/s
        r=(EN>1.5e-15).*(We.*N*2e-16.*exp(-7.248e-15./EN))+...
          (EN<=1.5e-15).*(We.*N*6.619e-17.*exp(-5.593e-15./EN));
    case 2 %two body attachment 1/s
       r=(EN>1.05e-15).*(We*N.*(8.889e-5*EN+2.567e-19))+...
         (EN<=1.05e-15).*(EN>=4.7512e-16).*(We*N.*(6.089e-4*EN-2.893e-19))+...
         (EN<4.7512e-16).*0;
    case 3 %three body attachment 1/s
       r=(EN>1.0e-19).*(We.*N*N*4.7778e-59.*EN.^-1.2749)+(EN<1.0e-19)*0.;                  
     case 4 %electron mobility m^2/V/s  
       r=(EN>=1.0e-19).*We*1e-2./e+(EN<1.0e-19)*1.539*N0/N; %this takes unphysical behaviour 
       %of Morrows mobility at low fields < 100 V/m at ground level
     case 5 %diffusion m^2/s 
        r=(EN>1.0e-19).*(0.3341e9*EN.^0.54069.*We./e*1e-2)+(EN<1.0e-19)*0.;     
     case 6 % Babaeva two body 1/s
         r=We.*N*4.3e-19.*exp(-1.05.*abs(5.3-log(EN/1e-17)).^3);        
      case 7 % Babaeva three body 1/s
         r=We.*N*N*1.6e-37.*(EN/1e-17).^-1.1 ;
        
    otherwise
        disp('Wrong input i')
end















