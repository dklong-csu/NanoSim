clc
clear variables
close all

%%
%   If you get an error that says you aren't finding the right GLIBCXX the
%   run the following in the terminal before launching Matlab:
%       export LD_PRELOAD=/lib/x86_64-linux-gnu/libstdc++.so.6

%%
settings.executable = "./GoldParticleSimulation";
settings.inp_file = "params.txt";
settings.diameter_file = "diams.txt";
settings.density_file = "density.txt";
settings.verbose = true;

model_prms = [9.999997e+07
3.881126e-01
3.956051e-01];

[d,q] = simulatePSDs(model_prms,settings);

figure
hold on
for ii=1:size(q,2)
    plot(d,q(:,ii),'-o')
end
hold off