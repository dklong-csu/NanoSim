clear
clc
close all

% psd_dense = readmatrix("step05-densesolve-2500-particles-PSD.txt");
% times_dense = readmatrix("step05-densesolve-2500-particles-solvetimes.txt");
% 
% figure
% plot(0.3 * (3:2500) .^(1/3), psd_dense)
% 
% figure
% histogram(times_dense)

size = 18000;
fname = "step05-sparsesolve-" + size + "-particles";
psd_sparse = readmatrix(fname + "-PSD.txt");
times_sparse = readmatrix(fname + "-solvetimes.txt");

figure
plot(0.3 * (3:size) .^(1/3), psd_sparse)
figure
plot(3:size, psd_sparse)
% figure
% solve_types = cell(1);
% solve_types{1} = 'sparse';
% data = cell(1);
% data{1} = times_sparse;
% violinplot(data,solve_types,'ShowBox',false);