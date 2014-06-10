%% Simulation SFCW Radar RCS estimation
% Bachelor thesis by Sebastian M�ller

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

R = 10;     % distance radar - target [m]
v = 0;      % speed of target [m/s] (v>0 -> moves towards receiver)
sigma = 5;  % radar cross section [m^2]

% radar constants

f_a = 1e6; % sampling frequency [Hz]

f_0 = 24.125e9;   % center frequency [Hz]
B = 50e6;     % sweep frequency [Hz]
T_f = 10e-3;  % sweep time per flank [s]
n = 64;      % frequency steps per flank [1]
N = 1;      % measuring intervals [1]

P_s = 100;  % transmission power [W]
G_T = 100;  % transmitting antenna gain [1]
G_R = 100;  % receiving antenna gain [1]

N_fft = 4096*128; % FFT Size

%% Calculate remaining data

T = N*2*T_f; % calculate measure time
t = 0:1/f_a:T-1/f_a; % create time vector
lambda = c / f_0; % wavelength

f_sfcw = SFCW_Freq(t, 2*T_f, f_0, B, n, N); % create sfcw frequency vector

tau = 2*R/c;    % time delay
f_D = 2*v*f_0/c;    % doppler frequency

f_min = f_0 - B/(2*n);

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

%% Diagnostic debug data

df = B/n;
dt = T_f/n;
R_max = c/(2*df)
dR = c/(2*B)

phi_bb = baseband_phase(f_min, tau, dt, df, T/N, t); % get phase of baseband signal

% a = angle(q);
% dp = a(:,100);
% dp = abs(dp - a(:,105+floor(dt*f_a)));
% R_est = dp*c/(4*pi*df);
%% simulation data

% transmitted signal

s = sqrt(P_s) * exp(1i * 2 * pi * discrete_int(f_sfcw, f_a, 0, T));

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

% q = c_r * P_s * exp(1i*2*pi*(discrete_int(f_sfcw, f_a, 0, T)...
%     - [zeros(1, timeToInt(tau, f_a)) discrete_int(f_sfcw, f_a, 0, T-tau)]...
%     - f_D * t));
t_dec = decimate(t, 300);
t_dec  = interp(t_dec, 300, 1);
t_dec = t_dec(1:length(t));
q = c_r * P_s * exp(1i*2*pi*(phi_bb - f_D * t_dec));
    % Filled time delayed vector with truncated values at beginning to
    % simulate periodizity

    
% add noise  
%q = awgn(q, 60);

%% Estimation

[q_fft, f_D_est, f_R_est] = spektrum(q, N_fft, f_a); % Calculate fft and frequency shift

% Calculate estimated v

v_est = c/2 * -f_D_est/f_0; % use negative f_D because we are measuring '-f_D' in analytical term

% Calculate estimated R

kappa = f_R_est; % calculate spatial frequency (still buggy)
R_est = c/(2*B)*kappa; % estimate R



%% Plot

%plot abs
fig = figure(1);
subplot(3,1,1);
plot(t, abs(q));
title('baseband signal abs');
xlabel('t/s');
ylabel('Ampl.');

% plot phase
subplot(3,1,2);
plot(t, angle(q));
title('baseband signal phase');
xlabel('t/s');
ylabel('Phase/rad');
subplot(3,1,3);

% FFT Plot
x_fa = 0:f_a/N_fft:f_a-f_a/N_fft;
plot(x_fa-f_a/2, q_fft, 'b.-')
%axis([-fn fn 0 (max(y)-min(y))/4*1.1])
title('FFT')
ylabel('Amplitude')
xlabel(['Aufl�sung: ',num2str(df),' Hz Frequenz in Hz'])
grid

% SFCW frequency plot
% figure(fig+1);
% plot(t, f_sfcw);
% title('SFCW Frequency');
% xlabel('t/s');
% ylabel('f/Hz');

% RV Plot
% figure(fig+1);
% R_vec = linspace(R-100, R+100, 1000); % generate R vector for RV plot
% v_vec = R_vec .* B/(f_0 * T_f) - c*kappa/(2*f_0*T_f); % generate v vector for RV plot
% plot(R_vec, v_vec, R_est, v_est, 'rx', R, v, 'gx');
% title('RV-Plot');
% xlabel('R / m');
% ylabel('v / m/s');
% legend('RV-Plot', 'Estimated', 'Actual');

%% print results
fprintf('Symbol \t\tValue\n');
fprintf(['v\t\t', num2str(v), ' m/s\n']);
fprintf(['v (estimated)\t', num2str(v_est), ' m/s\n']);
fprintf(['R\t\t', num2str(R), ' m\n']);
fprintf(['R (estimated)\t', num2str(R_est), ' m\n']);


