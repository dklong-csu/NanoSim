clc
clear variables
close all

%%
settings.executable = "./GoldParticleSimulation";
settings.inp_file = "params.txt";
settings.diameter_file = "diams.txt";
settings.density_file = "density.txt";
settings.verbose = false;

[data_diams, data_psds] = cleanRawData("Au_quench_corediameter_hplc.txt",1);

%%  Optimization settings
num_iters = 100;

F = @(x) costFcn(x,settings, data_diams, data_psds);

lb = [1,0,0];
ub = [1e10, 1,1];

x0 = [1e8,0.5,0.1];

%%  Global Optimization
globalopts = optimoptions('simulannealbnd',...
    'MaxIterations',num_iters,...
    'Display','iter');
xglobal = simulannealbnd(F,x0,lb,ub,globalopts);

%%  Local Optimization
localopts= optimoptions('fmincon',...
    'Display','iter-detailed',...
    'MaxIterations',num_iters);

xbest = fmincon(F,xglobal,[],[],[],[],lb,ub,[],localopts);

%%  Plot optimal
[dsim,qsim] = simulatePSDs(xbest,settings);
for ii=1:size(data_psds,2)
    figure
    hold on
    plot(data_diams, data_psds(:,ii));
    scatter(dsim, qsim(:,ii));
    xlim([0, 20])
    hold off
end
%%
function L = costFcn(prm,settings,data_diam,data_psds)
    [d,q] = simulatePSDs(prm,settings);
    dlo = d(1) - 1e-6;
    dhi = d(end) + 1e-6;
    dinterp = [dlo;d;dhi];
    L = 0;
    for ii=1:size(q,2)
        psd = [0;q(:,ii);0];
        interp = griddedInterpolant(dinterp,psd,'makima','nearest');
        sim = interp(data_diam);
        L = L + norm(sim - data_psds(:,ii));
    end
end