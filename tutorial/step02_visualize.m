clear
clc
close all

conc = readmatrix('step02-01-output-concentrations.txt');
t    = readmatrix('step02-01-output-times.txt');

figure
hold on
plot(t, conc(:,1), 'DisplayName',"[A]")
plot(t, conc(:,2), 'DisplayName',"[B]")
plot(t, conc(:,3), 'DisplayName',"[C]")
plot(t, conc(:,4), 'DisplayName',"[D]")
legend
xlabel('Time')
ylabel('Concentration')
title("First set of parameters")
hold off

conc = readmatrix('step02-02-output-concentrations.txt');
t    = readmatrix('step02-02-output-times.txt');

figure
hold on
plot(t, conc(:,1), 'DisplayName',"[A]")
plot(t, conc(:,2), 'DisplayName',"[B]")
plot(t, conc(:,3), 'DisplayName',"[C]")
plot(t, conc(:,4), 'DisplayName',"[D]")
legend
xlabel('Time')
ylabel('Concentration')
title("Second set of parameters")
hold off