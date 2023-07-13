clear
clc
close all

data = readmatrix('step03-data-conc.txt');
bad  = readmatrix('step03-badsim-conc.txt');
good = readmatrix('step03-goodsim-conc.txt');
t    = readmatrix('step03-times.txt');

figure
hold on
plot(t, bad(:,1),"Color","#0072BD")
plot(t, bad(:,2),"Color","#D95319")
plot(t, bad(:,3),"Color","#EDB120")
plot(t, bad(:,4), "Color","#7E2F8E")

scatter(t, data(:,1),"MarkerEdgeColor","#0072BD")
scatter(t, data(:,2),"MarkerEdgeColor","#D95319")
scatter(t, data(:,3),"MarkerEdgeColor","#EDB120")
scatter(t, data(:,4), "MarkerEdgeColor","#7E2F8E")
xlabel('Time')
ylabel('Concentration')
title("Bad parameters")
hold off

figure
hold on
plot(t, good(:,1),"Color","#0072BD")
plot(t, good(:,2),"Color","#D95319")
plot(t, good(:,3),"Color","#EDB120")
plot(t, good(:,4), "Color","#7E2F8E")

scatter(t, data(:,1),"MarkerEdgeColor","#0072BD")
scatter(t, data(:,2),"MarkerEdgeColor","#D95319")
scatter(t, data(:,3),"MarkerEdgeColor","#EDB120")
scatter(t, data(:,4), "MarkerEdgeColor","#7E2F8E")
xlabel('Time')
ylabel('Concentration')
title("Good parameters")
hold off