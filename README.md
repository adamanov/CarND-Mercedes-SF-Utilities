# SFND Radar Target Generation and Detection
## [Rubric](https://review.udacity.com/#!/rubrics/2548/view) Points
---
#### 1. FMCW Waveform Design
Using the given system requirements, design a FMCW waveform. Find its Bandwidth (B), chirp time (Tchirp) and slope of the chirp.

```Matlab
%% Radar Specifications 
% Frequency of operation = 77GHz
% Max Range = 200m
% Range Resolution = 1 m
% Max Velocity = 100 m/s
R_Target = 110;    % initial distance of the target 
v_Traget = -20;    % speed of the target. The velocity varies from -70 to +70 m/s

%Operating carrier frequency of Radar 
fc= 77e9;            % carrier freq            [Hz= 1/s]

% Initial conditions of Radar
max_range = 200;     % the radar maximum range [m]
d_res = 1;           % range resolution        [m]
c = 3*10^8 ;         % ligth speed             [m/s]

% Sweep time for each chirp is defined as rule by 5.5 times of
% round trip time for Maximum Range
Tchirp = 5.5 * 2 * max_range / c;       % Tchirp

% Bandwidth for the each chirp for given resolution 
Bandwidth = c/(2*d_res);

% The slope of the chirp
Slope = Bandwidth/Tchirp;
```
#### 2. Simulation Loop
Simulate Target movement and calculate the beat or mixed signal for every timestamp.

```Matlab
%The number of chirps in one sequence. Its ideal to have 2^ value for the ease of running the FFT
%for Doppler Estimation. 

Nd=128;                   % #of doppler cells OR #of sent periods % number of chirps
Nr=1024;                  %for length of time OR # of range cells
t=linspace(0,Nd*Tchirp,Nr*Nd); %total time for samples

%Creating the vectors for Tx, Rx and Mix based on the total samples input.
Tx  = zeros(1,length(t)); %transmitted signal
Rx  = zeros(1,length(t)); %received signal
Mix = zeros(1,length(t)); %beat signal

%Similar vectors for range_covered and time delay.
r_t = zeros(1,length(t));
td  = zeros(1,length(t));
for i=1:length(t)          
    % *%TODO* :
    %For each time stamp update the Range of the Target for constant velocity. 
    r_t(i) = R_Target + v_Traget*t(i);
    tau(i) = (2*r_t(i))/c;
    % *%TODO* :
    %For each time sample we need update the transmitted and
    %received signal. 
    Tx(i)  = cos(2*pi * (fc*t(i) +Slope*(t(i)^2) / 2));
    Rx(i) = cos(2*pi * (fc*(t(i)-tau(i)) + 0.5*Slope*(t(i)-tau(i))^2));
    
    % *%TODO* :
    %Now by mixing the Transmit and Receive generate the beat signal
    %This is done by element wise matrix multiplication of Transmit and
    %Receiver Signal
    Mix(i) = Tx(i).*Rx(i);
    
end
```

