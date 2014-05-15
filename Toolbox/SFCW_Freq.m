function [ f ] = SFCW_Freq( t, T, f_0, B, n, N )
%SFCW_Freq generates frequency vector for SFCW Radar
%   Input:
%   t:      time vector
%   T:      periodic time
%   f_0:    center frequency
%   B:      sweep frequency (symm. around f_0)
%   n:      steps per flank
%   N:      measuring intervals
%   Output:
%   f:      frequency vector
%
%   Example: SMCW_Freq()

dt = T/(2*n); % time step
df = B/n; % frequency step
f = f_0 - B/2; % start frequency (minimal freq)

% measuring intervals loop
for j=0:1:N-1
    %% rising flank
    
    % step loop
    for i=j*2*n+1:1:(2*j+1)*n
        f = f + df*heaviside(t-i*dt);
    end

    %% falling flank

    %step loop
    for i=(2*j+1)*n+1:1:(j+1)*2*n
       f = f - df*heaviside(t-i*dt); 
    end
end
%plot(t, f)

end

