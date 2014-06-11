function [y_fft, fR_x] = spektrum(y, N, f_a) 
%spektrum generates fft of signal and finds abs of highest peak
%   Input:
%       y:      signal vector
%       N:      N-point-fft
%       f_a:    sampling frequency
%   Output:
%       y_fft:  fft vector of signal
%       fR_x:   abs of highest peak

% do fft

df = f_a/N;
fn = f_a/2;
x_fa = 0 : df : f_a-df;

H = fft(y, N);
amplH = abs(H);
amplitudengang = fftshift(amplH/N);

% find peaks

[peak_y, peak_x] = findpeaks(amplitudengang, 'NPEAKS', 1, 'SORTSTR', 'descend'); ...
	% finds highest peak in spectrum

fR_x = abs(x_fa(:,round(peak_x))-fn); % Get Highest peak freqency

y_fft = amplitudengang; % assign fft vector

end 