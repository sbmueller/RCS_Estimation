function [ f ] = SMCW_Freq( t, T, f_s0, B, N, K )
%SMCW_Freq_alt Generiert Frequenz- und Zeitvektor für SMCW-Radar
%   Eingabe:
%   t:      Zeitvektor
%   T:      Länge einer Periode
%   f_s0:   Mittenfrequenz
%   B:      Frequenzhub (symm. um f_s0)
%   N:      Stepanzahl pro Flanke, sollte gerade sein
%   K:      Periodenanzahl
%   Ausgabe:
%   f: Frequenzvektor
%
%   Beispiel: SMCW_Freq()

dt = T/(2*N); % Dauer eines Steps
df = B/N; % Bandbreite eines Steps
f = f_s0 - B/2; % Startfrequenz - Minimale Frequenz

for j=0:1:K-1
    %% Steigende Flanke
    for i=j*2*N+1:1:(2*j+1)*N
        f = f + df*heaviside(t-i*dt);
    end

    %% Fallende Flanke

    for i=(2*j+1)*N+1:1:(j+1)*2*N
       f = f - df*heaviside(t-i*dt); 
    end
end
%plot(t, f)

end

