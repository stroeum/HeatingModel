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
flag.plot         = 'on';
flag.movie        = 'off';

U          = 19.2;
h          = 0;
tau0Dn0km  = Heating1D( U, h, flag); % 1-D
hgexport(1,'1D-0km');