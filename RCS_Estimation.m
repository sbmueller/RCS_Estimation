%% Simulation FMCW Radar
% Bachelorarbeit von Sebastian Müller

% Lösche alle Variablen
clear all

% Lösche Konsole
clc;

disp('// SIMULATION FMCW RADAR')
disp(' ')
addpath Toolbox/
%% Parameter
% Naturkonstanten

c = 299792458; % Lichtgeschwindigkeit [m/s]
t = linspace(0, 1, 10000);

% Objektkonstanten

R = 10; % Distanz zum Ziel [m]
v = 0; % Geschwindigkeit Ziel [m/s] (v>0 -> Auf Empfänger zu)
sigma = 5; % Radarquerschnitt in [m^2]

% Radardaten

f_s0 = 10; % Mittenfrequez [Hz]
B = 20; % Frequenzhub Bandbreite [Hz]
T_f = 0.5; % Flankendauer [s]
n = 5; % Stepanzahl pro Flanke [1]
N = 1; % Messintervalle [1]


% Erstellen der übrigen Daten

f_smcw = SMCW_Freq(t, 2*T_f, f_s0, B, n, N);
f = @(time) f_smcw(:, t==time);
%plot(t, f(t))

%% Simulationsdaten

%Sendesignal
s = cos(2*pi*f(t).*t);
ft = abs(fft(s));
plot (ft);