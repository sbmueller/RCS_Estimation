function [ int ] = discrete_int( y, dx, min, max )
%discrete_int integrate discrete signals by adding up areas from 0 to max
%   Input:
%   y:      function vector
%   dx:     step interval (smallest time resolution)
%   max:    upper bound
%   Output:
%   int:    integrated vector

int = zeros(1, floor(max)-floor(min));
temp = 0;
for i=floor(min):1:floor(max-1)
    temp = temp + y(:,i)*dx;
    int(:,i-floor(min)+1) = temp;
end

end

