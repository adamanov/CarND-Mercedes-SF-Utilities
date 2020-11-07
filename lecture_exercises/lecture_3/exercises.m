%  Content: 
%
%  1- Range   Estimation
%  2- Doppler Estimation
%  3- 1D-FFT 
%  4- 2D-FFT

%% Range estimation equation 
clear ; close all ; clc;
% R = c*Ts*fb / ( 2*Bsweep)  

%Ts = chirp Time, Bsweep = chirp Bandwidth
% Those values are determined as we define the configuration of radar
% based on its range resolution and trip time for Radars maximum range



% Beat frequency 
% fb = f_ramping - f_received

% -----------------------------------------------------%

% The range resolution = 1m 
%

range_max = 300;    % The radar maximum range = 300m 
c = 3*10^8 ;        % ligth speed [m/s]
d_res = 1 ;         % range resolution [m]  d_res = c/(2*B_sweep) "higher Bandwidth the smaller range resoulation" 


% signal trip time t= 2*range/c


% TODOs to calculate the range in meters of four targets 
% with respective measured beat frequencies [0 MHz, 1.1 MHz, 13 MHz, 24 MHz]. 

% TODO : Find the Bsweep of chirp for 1 m resolution

B_sweep = c/2 * d_res;

% TODO : Calculate the chirp time based on the Radars Max Range
% T_chirp = 5.5 * 2 *range_max /c

T_s = 5.5 * (2 * range_max) /c ;

% TODO: define the frequncy shift 
fb = [ 0 , 1.1e6 13e6 24e6];
calcualted_range = c * T_s *fb / (2*B_sweep);  %[m]

% display the calcualte range

disp(calcualted_range);


%% Doppler Estimation 
clear; close; clc;
c = 3*10^8 ; % the speed of light

f = 77*10^9; % [Hz] the daras operating freq.

% TODO: Calculate the wavelength

lambda = c/f ;

% TODO: Define the doppler shifts in Hz using 

doppler_shift = [3e3 4.5e3 11e3 -3e3];

% Calculate the velocity of the targets fd= 2*vr/lambda

vr = doppler_shift*lambda /2 ;

disp(vr)


%% Fast Fourier Transform Exercise 
clear; close; clc;

Fs = 1000;       % Sampling frequency 
T = 1/Fs;        % Sampling period
L = 1500;        % Length of signal
t = (0:L-1)*T;   % Time vector 

% Form a signal containing a 77Hz sinusoid of amplitude 2.
A = 2;
f = 77;

S = A*cos(2*pi*f*t);
% Form a signal containing a 77 Hz sinusoid of ampilutude 0.7 and
% a 43 Hz sinusoid of amplitude 2.
% S = 0.7*sin(2*pi*77*t) + 2*sin(2*pi*43*t);

% Corrpt the signal with noise 

X = S + 2*rand(size(t));

% Plot the noisy signal in the time domain. 
% It is diffucult to identify the frequncy components 
% by looking at signal X(t);
plot (1000*t(1:50) , X(1:50));
title('Signal COrrupted with zero-mean random noise')
xlabel('t (milliseconds)')
ylabel('X(t)')

% TODO: Compute the Fourier transform of the signal

signal_fft = fft(X);

% TODO: Compute the two sided spectrum P2. 
% Then compute the single-sided spectrum P1 based on P2 
% and the even- valued signal length L.

P2 = abs(signal_fft/L); % take the amplitude of the normalized signal
P1 = P2(1:L/2+1); % 

% Ploting 
f = Fs * (0:(L/2))/L;
plot(f,P1);
title('Single -Sided Amplitude Spectrum of X(t)')
xlabel('f [Hz]')
ylabel('[P1(f)]')

%%  2D FFT Exercise
clear; close; clc;

% 2-D transform 
% The 2-D Fourier transform is useful for processing 2-D signals and other
% 2-D data such as images.
% Create and plot 2-D data with repeated blocks

P = peaks(20); 
X = repmat(P,[5 10]);
imagesc(X);

% Compute the 2-D Fourier transform of the data
% Shift the zero-frequency component to the center of the output, 
% and plot the resulting 100-by-200 matrix, which is the same size as X.

Y = fft2(X);
imagesc(abs(fftshift(Y)))

% Pad X with zeros to compute a 128-by-256 transform.
Y = fft2(X,2^nextpow2(100),2^nextpow2(200));
imagesc(abs(fftshift(Y)));



