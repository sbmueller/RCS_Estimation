function [ int ] = discrete_int( y, f_a, min, max )
%discrete_int integrate discrete signals by adding up areas from 0 to max
%   Input:
%   y:      function vector
%   dx:     step interval (smallest time resolution)
%   max:    upper bound
%   Output:
%   int:    integrated vector

%% Using for loop (previous version, does not work anymore)
% int = zeros(1, floor(max)-floor(min));
% temp = 0;
% for i=floor(min):1:floor(max-1)
%     temp = temp + y(:,i)*dx;
%     int(:,i-floor(min)+1) = temp;
% end

%% Using cumtrapz() function (better)

% convert time parameters to index parameters
min = floor((f_a * min)+1);
max = floor(f_a * max);

int = cumtrapz(y(min:max))/f_a; % see cumtrapz documentation

end

