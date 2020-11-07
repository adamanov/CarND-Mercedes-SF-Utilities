clc, clear, close all

% Operating frequency [Hz] or [1/s]
fc = 77.0e9;

% Transmitted Power [W]
Pt = 3e-3;

% Antenna Gain [linear] or [dBi]
G = 10000;

% Minimum  Detectable Power
Ps = 1e-10;

% RCS of a car [m^2]
RCS = 100;

% Speed of light  [m/s]
c = 3*10^8;

% TODO: Calculate the wavelength [m]

lambda = c/fc;

% TODO: Measure the Maximum Range a Radar can see.
%
% Range equation -> R
% =(Pt*power(G,2)*power(lambda,2)*power(RCS,2)/(Ps*power(4/pi,3))

nom = Pt*power(G,2)*power(lambda,2)*RCS;
den = Ps*power(4*pi,3);

range= power(nom/den,1/4);
disp(range);



