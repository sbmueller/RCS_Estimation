function [y_fft, f_fft] = spektrum(y, N, f_a) 
%spektrum generates the spectrum and frequency vector of signal y
%   Input:
%       y:      signal vector
%       N:      N-point-fft
%       f_a:    sampling frequency
%   Output:
%       y_fft:  fft vector of signal
%       f_fft:	frequency vector

% do fft

df = f_a/N;
fn = f_a/2;
x_fa = 0 : df : f_a-df;

H = fft(y, N);
amplH = abs(H);
amplitudengang = fftshift(amplH/N);

y_fft = amplitudengang; % assign fft vector
f_fft = x_fa-fn;

% plot(f_fft, y_fft);
end 