function spektrum(y, N, f_a) 
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

stem(x_fa-fn, amplitudengang, 'b.-')
%axis([-fn fn 0 (max(y)-min(y))/4*1.1])
title('Amplitudengang')
ylabel('Amplitude')
xlabel(['Auflösung: ',num2str(df),' Hz Frequenz in Hz'])
grid
end 