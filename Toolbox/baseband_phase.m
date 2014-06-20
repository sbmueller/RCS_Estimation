function [ f_int ] = baseband_phase( f_0, tau, dt, df, T )
%baseband_phase creates the phase of the baseband signal with the
%analytical term
%   f_0:    minimum frequency
%   tau:    radar time delay
%   dt:     sweep time
%   df:     sweep frequency
%   T:      interval time


t_jump = dt:dt:T;    % create vector of jump times

f_intup = tau*(f_0 + t_jump(1:floor(length(t_jump)/2 + 1)) * df/dt); % Integral of rising flank

f_intdown = tau * (-t_jump(floor(length(t_jump)/2 + 2):length(t_jump)) * df/dt) ...
	+ 2*f_intup(:, length(f_intup)) - f_0*tau; % Integral of falling flank

f_int = 2*pi* [f_intup f_intdown]; % combine both flanks

%% Upsampling

% rep_factor = floor(length(t)/length(f_int));
% % interpolate phase to fit time vector with zero order function
% f_int = repmat(f_int, rep_factor, 1);
% f_int = reshape(f_int, 1, rep_factor*length(t_jump));

% % repeat last few samples to fit odd number of timesamples
% for i=1:length(t)-length(f_int)
%     f_int = [f_int f_int(:,length(f_int))];
% end

% debug plot
% stairs(t_jump, f_int);
