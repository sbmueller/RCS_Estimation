clear;
t = linspace(0,1,1000);
q = exp(1i * 2*pi* 100 .*t);

q_fft = fft(q);

plot(abs(fft(q)));