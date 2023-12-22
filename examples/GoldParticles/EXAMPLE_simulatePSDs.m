clc
clear variables
close all

%%
settings.executable = "./GoldParticleSimulation";
settings.inp_file = "params.txt";
settings.diameter_file = "diams.txt";
settings.density_file = "density.txt";
settings.verbose = true;

model_prms = [8.54958e+06
0.5
0.01];

[d,q] = simulatePSDs(model_prms,settings);

figure
hold on
for ii=1:size(q,2)
    plot(d,q(:,ii),'-o')
end
hold off