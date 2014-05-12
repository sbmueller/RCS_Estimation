function [ y ] = heaviside( x )
%heaviside Implementation des Einheitssprungs
%   Definiert als   0 für x < 0
%                   1 für x >= 0
y = (x >= 0);


end

