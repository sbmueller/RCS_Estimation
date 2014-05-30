function [ i ] = timeToInt( t, f_a )
%timeToInt convert time value to corresponding index value in time vector
%   Input:
%   t:      time value
%   f_a:    sampling frequency
%   Output:
%   i:      index value

    i = round(f_a * t);