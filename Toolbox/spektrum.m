function [y_fft, f_peak] = spektrum(y, N, f_a) 
%spektrum generates spectral plot of signal y
%   Input:
%       y:      signal vector
%       N:      N-point-fft
%       f_a:    sampling frequency

df = f_a/N;
fn = f_a/2;
x_fa = 0 : df : f_a-df;

H = fft(y, N);
amplH = abs(H);
amplitudengang = fftshift(amplH/N);

[peak_y, peak_x] = findpeaks(amplitudengang, 'MINPEAKHEIGHT', ...
    0.3*(max(amplitudengang) - min(amplitudengang)));
peak_x = round(mean(peak_x));
f_peak = x_fa(:,peak_x)-fn;
y_fft = amplitudengang;

plot(x_fa-fn, amplitudengang, 'b.-')
%axis([-fn fn 0 (max(y)-min(y))/4*1.1])
title('FFT')
ylabel('Amplitude')
xlabel(['Auflösung: ',num2str(df),' Hz Frequenz in Hz'])
grid
end 