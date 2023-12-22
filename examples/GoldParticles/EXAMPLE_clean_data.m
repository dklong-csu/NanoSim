clc
clear variables 
close all

%%
[diams, PSDs] = cleanRawData("Au_quench_corediameter_hplc.txt",1);

figure
plot(diams, PSDs)
xlim([0 20])