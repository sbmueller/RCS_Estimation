function [y_fft, fD_x, fR_x] = spektrum(y, N, f_a) 
%spektrum generates fft of signal and finds abs of highest peak
%   Input:
%       y:      signal vector
%       N:      N-point-fft
%       f_a:    sampling frequency
%   Output:
%       y_fft:  fft vector of signal
%       fD_x:   (will be removed) Doppler peak
%       fR_x:   abs of highest peak

df = f_a/N;
fn = f_a/2;
x_fa = 0 : df : f_a-df;

H = fft(y, N);
amplH = abs(H);
amplitudengang = fftshift(amplH/N);

[peak_y, peak_x] = findpeaks(amplitudengang, 'MINPEAKHEIGHT', ...
    0.5*(max(amplitudengang) - min(amplitudengang)), 'SORTSTR', 'descend', ...
    'NPEAKS', 6); % Find two highest peaks in fft, f_D is in the middle

fD_x = round(mean(peak_x)); % Find middle of two peaks
fD_x = x_fa(:,fD_x)-fn; % Get corresponding frequency of peak

[peak_y, peak_x] = findpeaks(amplitudengang, 'NPEAKS', 1, 'SORTSTR', 'descend');
fR_x = abs(x_fa(:,round(peak_x))-fn); % Get Highest peak freq

y_fft = amplitudengang;

% plot


end 