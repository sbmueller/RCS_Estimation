%% Simulation SFCW Radar RCS estimation
% Bachelor thesis by Sebastian Müller

% clear all variables
clear all

% clear console
clc;

disp('**** SIMULATION SFCW RADAR ****')
disp(' ')
addpath Toolbox/

%% Parameters
% natural constants

c = 299792458; % speed of light [m/s]

% object constants

R = 100;     % distance radar - target [m]
v = 0;      % speed of target [m/s] (v>0 -> moves towards receiver)
sigma = 5;  % radar cross section [m^2]

% radar constants

f_a = 100e3; % sampling frequency [Hz]

f_0 = 24.125e9;   % center frequency [Hz]
B = 15e6;     % sweep frequency [Hz]
T_f = 10e-3;  % sweep time per flank [s]
n = 128;      % frequency steps per flank [1]
N = 1;      % measuring intervals [1]

P_s = 100;  % transmission power [W]
G_T = 100;  % transmitting antenna gain [1]
G_R = 100;  % receiving antenna gain [1]

%% Calculate remaining data

T = N*2*T_f; % calculate measure time
S = floor(f_a*T);   % number of samples
t = 0:1/f_a:T-1/f_a; % create time vector
lambda = c / f_0; %! CHANGES IN TIME

f_sfcw = SFCW_Freq(t, 2*T_f, f_0, B, n, N); % create frequency vector

tau = 2*R/c;    % time delay
f_D = 2*v*f_0/c;    % doppler frequency

%% check values

if R > c/2 * n/B
    disp(['WARNING: R too big for unambiguous calculation. R has to be smaller than ' , num2str(c/2*n/B) , ' m']);
    disp(' ');
end
if R < c/2 * 1/B
    disp('WARNING: R too small for radial resolution. Use bigger sweep frequency or increase R.');
    disp(' ');
end
if f_a/2 < B/n
    disp('WARNING: Sampling freq too small for sweep freq!');
    disp(' ');
end
if f_a/2 < f_D
    disp('WARNING: Sampling freq too small for doppler freq!');
    disp(' ');
end
if v < c/(2*f_0*T_f)
    disp('WARNING: v too small for speed resolution! Use longer sweep time.');
    disp(' ');
end
%% simulation data

% transmitted signal

s = sqrt(P_s) * exp(1i * 2 * pi * discrete_int(f_sfcw, 1/f_a, 1, S));

c_r = (sqrt(G_T)*sqrt(G_R)*lambda*sqrt(sigma))/((4*pi)^(3/2)*R^2); ...
% attenuation because of radar equation

% Check possibility
if c_r > 1
    c_r = 1;
end

% received signal

%% Simulate mixing (not recommended, not further developed)
% r = c_r * sqrt(P_s) * exp(1i * 2 * pi *...
%     discrete_int(f_sfcw, 1/f_a, 1, floor(f_a*(T-tau))));
% r = [zeros(1, S - floor(f_a*(T-tau))) r]; % fill received vector with nulls for same length
% r = r .* exp(1i * 2 * pi * f_D * t); % add phase shift because of doppler
% 
% % mixing
% 
% q = s.*conj(r);

%% OR use analytical term (recommended)

q = c_r * P_s * exp(1i*2*pi*(discrete_int(f_sfcw - f_D, 1/f_a, 1, S)...
    - [discrete_int(f_sfcw, 1/f_a, f_a*(T-tau), S) discrete_int(f_sfcw, 1/f_a, 1, f_a*(T-tau))]));...
    % Filled time delayed vector with truncated values at the beginning to
    % simulate periodicity
% 
% q = c_r * P_s * exp(1i*2*pi*(discrete_int(f_sfcw, tau, floor(f_a*(T-tau)), S)...
%     - f_D*t));
    
% add noise  
%q = awgn(q, 60);
%% Plot

% plot abs
fig = figure(1);
subplot(3,1,1);
plot(t(1:S-1), abs(q));
title('baseband signal abs');
xlabel('t/s');
ylabel('Ampl.');

% plot phase
subplot(3,1,2);
plot(t(1:S-1), unwrap(angle(q)));
title('baseband signal phase');
xlabel('t/s');
ylabel('Phase/rad');
subplot(3,1,3);

% plot fft
[q_fft, f_D_est] = spektrum(q, 4096, f_a);

% SFCW frequency plot
% figure(fig+1);
% plot(t, f_sfcw);
% title('SFCW Frequency');
% xlabel('t/s');
% ylabel('f/Hz');

%% Estimation

% Calculate estimated v

v_est = c/2 * -f_D_est/f_0; % use negative f_D because we are measuring '-f_D' in analytical term

%% print results
fprintf('Symbol \t\tValue\n');
fprintf(['v\t\t', num2str(v), ' m/s\n']);
fprintf(['v (estimated)\t', num2str(v_est), ' m/s\n']);
fprintf(['R\t\t', num2str(R), ' m\n']);
a = angle(q);
a = a(1:length(a)-1) - a(2:length(a));
a = abs(a(:,8));
disp(' ');
disp(['Phasejump: ', num2str(a)]);
