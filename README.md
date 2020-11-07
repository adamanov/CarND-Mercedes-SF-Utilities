# SFND Radar Target Generation and Detection
## [Rubric](https://review.udacity.com/#!/rubrics/2548/view) Points
---
#### 1. FMCW Waveform Design
Using the given system requirements, design a FMCW waveform. Find its Bandwidth (B), chirp time (Tchirp) and slope of the chirp.
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
