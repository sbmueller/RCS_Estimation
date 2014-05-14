function [ y ] = heaviside( x )
%heaviside Implementation of heaviside step function
%   Defined as      0 for x < 0
%                   1 for x >= 0
y = (x >= 0);


end

