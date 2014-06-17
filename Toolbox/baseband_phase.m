function [ f_int ] = baseband_phase( f_0, tau, dt, df, T, t )
%baseband_phase creates the phase of the baseband signal with the
%analytical term
%   f_0:    minimum frequency
%   tau:    radar time delay
%   dt:     sweep time
%   df:     sweep frequency
%   T:      interval time
%   t:      time vector


t_jump = 0:dt:T-dt;    % create vector of jump times
f_intup = f_0*tau-0.5*(-0.5* t_jump(1:floor(length(t_jump)/2 + 1)) *tau+tau^2)*df/dt; % integral of raising flank

f_intdown = 0.5*(-0.5* t_jump(floor(length(t_jump)/2 + 2):length(t_jump)) *tau+tau^2)*df/dt ...
    + 2*f_intup(:,length(f_intup)) - f_0*tau;
   % integral of falling flank

f_int = [f_intup f_intdown]; % combine both flanks

rep_factor = floor(length(t)/length(f_int));
% interpolate phase to fit time vector with zero order function
f_int = repmat(f_int, rep_factor, 1);
f_int = reshape(f_int, 1, rep_factor*length(t_jump));

% repeat last few samples to fit odd number of timesamples
for i=1:length(t)-length(f_int)
    f_int = [f_int f_int(:,length(f_int))];
end

% debug plot
% plot(t, f_int);
