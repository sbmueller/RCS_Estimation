function [ sigma_est ] = estimate_sigma( q, P_s, param, R_est )
%estimate_sigma Core function of the RCS estimation
%   Estimates RCS via received power and estimated R
%   q:          baseband signal
%   P_s:        trasmit power
%   param:      constant eqals (4*pi)^3/(G_R*G_T*lambda^2)
%   R_est:      estimated R
%   sigma_est:  estimated RCS


% FFT is the best way to estimate amplitude of signal (kay)
q_fft =2* fft(q, 4096)/length(q);


[P_e, loc] = findpeaks(abs(q_fft), 'NPEAKS', 1, 'SORTSTR', 'descend'); % find peak of fft
    
q_fft = abs(q_fft(loc-4:loc+4)); % get 9 values next to peak

f_low = 0:8;
f_high = 0:9/256:8; % create interpolation vectors

q_fft = interp1(f_low, q_fft, f_high, 'PCHIP'); % interpolate fft from 9 to 256 samples

P_e = findpeaks(q_fft, 'NPEAKS', 1, 'SORTSTR', 'descend')^2/P_s; % Find peak wich equals amplitude of baseband signal sqrt(P_s)*sqrt(P_e)

sigma_est = P_e/P_s *param * R_est^4; % Use radar equation to estimate RCS

end