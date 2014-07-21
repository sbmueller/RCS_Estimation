function [ sigma_est ] = estimate_sigma( q, f_a, P_s, param, R_est )
%estimate_sigma Core function of the RCS estimation
%   Detailed explanation goes here

N = 4096;

% FFT is the best way to estimate amplitude of signal (kay)
q_fft =2* fft(q, N)/length(q);
P_e = findpeaks(abs(q_fft), 'NPEAKS', 1, 'SORTSTR', 'descend')^2/P_s;

sigma_est = P_e/P_s *param * R_est^4;

end

