%% Simulation SFCW Radar RCS estimation
% Bachelor thesis by Sebastian Müller

% clear all variables
clear all

% clear console
clc;

disp('// SIMULATION SFCW RADAR')
disp(' ')
addpath Toolbox/

%% Parameters
% natural constants

c = 299792458; % speed of light [m/s]

% object constants

R = 1000;     % distance radar - target [m]
v = 10;      % speed of target [m/s] (v>0 -> moves towards receiver)
sigma = 5;  % radar cross section [m^2]

% radar constants

f_a = 200e4; % sampling frequency [Hz]

f_0 = 2.4e9;   % center frequency [Hz]
B = 100e3;     % sweep frequency [Hz]
T_f = 1e-3;  % sweep time per flank [s]
n = 5;      % frequency steps per flank [1]
N = 1;      % measuring intervals [1]

P_s = 100;  % transmission power [W]

% calculate remaining data

T = N*2*T_f; % calculate measure time
S = uint64(f_a*T);   % number of samples
t = linspace(0, T, S); % create time vector

f_smcw = SMCW_Freq(t, 2*T_f, f_0, B, n, N); % create frequency vector

%% simulation data

tau = 2*R/c;    % time delay
f_D = 2*v*f_0/c;    % doppler frequency

% transmitted signal

s = sqrt(P_s) * exp(1i * 2 * pi * discrete_int(f_smcw, 1/f_a, S));

% received signal

r = 0.5 * sqrt(P_s) * exp(1i * 2 * pi * discrete_int(f_smcw, 1/f_a, ceil(f_a*(T-tau))));
r = [zeros(1, uint32(f_a*T - ceil(f_a*(T-tau)))) r]; % fill received vector with nulls for same length
r = r .* exp(1i * 2 * pi * f_D * t); % add phase shift because of doppler


plot (t, real(s), t, real(r));
legend('transmitted', 'received');
xlabel('t/s');
ylabel('Ampl.');



