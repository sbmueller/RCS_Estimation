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

R = 1000;     % distance radar - target [m]
v = 10;      % speed of target [m/s] (v>0 -> moves towards receiver)
sigma = 5;  % radar cross section [m^2]

% radar constants

f_a = 1e6; % sampling frequency [Hz]

f_0 = 24.125e9;   % center frequency [Hz]
B = 125e6;     % sweep frequency [Hz]
T_f = 10e-3;  % sweep time per flank [s]
n = 64;      % frequency steps per flank [1]
N = 4;      % measuring intervals [1]

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

%% simulation data

% transmitted signal

s = sqrt(P_s) * exp(1i * 2 * pi * discrete_int(f_sfcw, 1/f_a, 1, S));

c_r = (sqrt(G_T)*sqrt(G_R)*lambda*sqrt(sigma))/((4*pi)^(3/2)*R^2); % attenuation because of radar equation

% Check possibility
if c_r > 1
    c_r = 1;
end

% received signal

%% Simulate mixing (not recommended)
% r = c_r * sqrt(P_s) * exp(1i * 2 * pi *...
%     discrete_int(f_sfcw, 1/f_a, 1, floor(f_a*(T-tau))));
% r = [zeros(1, S - floor(f_a*(T-tau))) r]; % fill received vector with nulls for same length
% r = r .* exp(1i * 2 * pi * f_D * t); % add phase shift because of doppler
% 
% % mixing
% 
% q = s.*conj(r);

%% OR use analytical term (recommended)

f_delay = [f_sfcw(ceil(f_a*(T-tau)):length(f_sfcw)) f_sfcw(1:floor(f_a*(T-tau)))];    % frequency vector of received signal (delayed)

q = c_r * P_s * exp(1i*2*pi*(discrete_int(f_sfcw - f_delay - f_D, 1/f_a, 1, S)));  % filled time delayed terms truncated values at beginning

% plot

NFFT = 2^nextpow2(S);
Y = fft(q, NFFT)/S;
f = f_a/2*linspace(0,1,NFFT/2+1);

subplot(2,1,1);
%plot(f, 2*abs(Y(1:NFFT/2+1)));
plot(t, f_sfcw, t, f_delay);
legend('freq sent');
xlabel('Hz');
ylabel('Ampl.');
subplot(2,1,2);
plot(t, real(q));
legend('baseband signal');
xlabel('t/s');
ylabel('Ampl.');
