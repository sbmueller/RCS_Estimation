function [ f ] = SMCW_Freq_alt( T, f_s0, B, N )
%SMCW_Freq_alt Generiert Frequenz- und Zeitvektor für SMCW-Radar
%   Eingabe:
%   T:      Länge einer Periode
%   f_s0:   Mittenfrequenz
%   B:      Frequenzhub (symm. um f_s0)
%   N:      Stepanzahl, sollte gerade sein
%   Ausgabe:
%   f: Frequenzvektor für eine SMCW-Periode
%
%   Beispiel: SMCW_Freq()

N = N+1;
df = B/N;
t_helper = linspace(0, T, 3);
f_min = f_s0 - df*N/2;
f_max = f_s0 + df*N/2;

f_helper = [f_min f_max f_min];

t_lin = linspace(0,T,N);
f_lin = interp1(t_helper, f_helper, t_lin);

[t, f] = stairs(t_lin, f_lin);
%plot(t, f)

end

