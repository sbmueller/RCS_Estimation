clear
set(gca, 'FontSize', 14);
t = linspace(0,10,1000);
t0 = [t t+10 t+20];
f = 0.8.*t;
f0 = [f f f];

set(gca,'defaulttextinterpreter','latex');
plot(t0, f0)
xlabel('t');
ylabel('f');
set(gca, 'XTickLabel', ['0'; ' '; 'T'; ' '; '2T'; ' '; '3T'; ' '; '4T'; ' ';'5T']);
%set(gca, 'YTickLabel', []);
set(gca, 'YTickLabel', ['f_{min}'; ' '; ' '; ' '; ' '; ' '; ' '; ' '; 'f_{max}']);
