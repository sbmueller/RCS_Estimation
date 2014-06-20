function Hpsd = LDS(y, N, f_a) 
%LDS generates the power spectral density and frequency vector for a signal y
%   Input:
%       y:      signal vector
%       N:      N-point-fft
%       f_a:    sampling frequency
%   Output:
%       y_psd:  psd vector
%       f_psd:	frequency vector

Pxx = abs(fft(y, N)).^2/length(y)/f_a; 
Hpsd = dspdata.psd(Pxx(1:length(Pxx)/2), 'f_a', f_a);
end 