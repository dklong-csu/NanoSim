clear
clc
close all

conc = readmatrix('step01-output-concentrations.txt');
t    = readmatrix('step01-output-times.txt');

figure
hold on
plot(t, conc(:,1), 'DisplayName',"[A]")
plot(t, conc(:,2), 'DisplayName',"[B]")
plot(t, conc(:,3), 'DisplayName',"[C]")
plot(t, conc(:,4), 'DisplayName',"[D]")
legend
xlabel('Time')
ylabel('Concentration')
hold off