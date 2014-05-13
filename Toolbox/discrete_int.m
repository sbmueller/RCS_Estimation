function [ int ] = discrete_int( y, dx, max )
%discrete_int integrate discrete signals by adding up areas from 0 to max
%   Input:
%   y:      function vector
%   dx:     step interval (smallest time resolution)
%   max:    upper bound
%   Output:
%   int:    integrated vector

int = zeros(1, uint16(max));
temp = 0;
for i=1:1:max
    temp = temp + y(:,i)*dx;
    int(:,i) = temp;
end

end

