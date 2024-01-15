clc
clear variables
close all

%%  Start a diary to log the optimization results
diaryFileName = "optimizationLog.txt";
%   Append to previous version of log (true) or delete log first (false)
addToLog = false;
if ~addToLog
    %   Opening a file in write mode deletes its contents
    fileID = fopen(diaryFileName,'w');
    %   Now close it so we don't accidentally edit the file
    fclose(fileID);
end
diary(diaryFileName)

%%  Show how long each step takes
tstart = tic;

%%  Seed the rng for reproducibility
seedValue = 0;
rng(seedValue);
fprintf("Rng seeded with value: %f\n", seedValue);

%%
settings.executable = "./GoldParticleSimulation";
settings.inp_file = "params.txt";
settings.diameter_file = "diams.txt";
settings.density_file = "density.txt";
settings.verbose = false;

[data_diams, data_psds] = cleanRawData("Au_quench_corediameter_hplc.txt",1);

%%  Optimization settings
num_iters = 1000;

F = @(x) costFcn(x,settings, data_diams, data_psds);

lb = [1,0,0];
ub = [1e10, 1,1];

x0 = [1e8,0.5,0.1];

%%  Global Optimization
globalopts = optimoptions('simulannealbnd',...
    'MaxIterations',num_iters,...
    'Display','iter');
xglobal = simulannealbnd(F,x0,lb,ub,globalopts);

tendG = toc(tstart);
fprintf("\nGlobal optimization performed in: ")
prettyTimerDisplay(tendG)
fprintf("\n\n")

%%  Local Optimization
tstartL = tic;
localopts= optimoptions('fmincon',...
    'Display','iter-detailed',...
    'MaxIterations',num_iters);

xbest = fmincon(F,xglobal,[],[],[],[],lb,ub,[],localopts);

%%  Solve for optimal parameters
[dsim,qsim] = simulatePSDs(xbest,settings);

tendL = toc(tstartL);
fprintf("\nLocal optimization performed in: ")
prettyTimerDisplay(tendL)
fprintf("\n\n")

save("optimized","xglobal","xbest","dsim", "qsim",'-mat');
%%  Plot optimal
close all % temporary
%   Interpolate simulation
simInterp = griddedInterpolant(dsim, qsim(:,4),'makima');
simPlotPSD = max(0,simInterp(data_diams));
figure
hold on
box on
plot(data_diams, data_psds(:,4), ...
    'LineWidth',2,...
    'DisplayName','Data');
plot(data_diams, simPlotPSD,...
    'LineWidth',2,...
    'DisplayName','Simulation');
xlim([0,20])
set(gca,'linewidth',2)
legend
hold off

% for ii=1:size(data_psds,2)
%     figure
%     hold on
%     plot(data_diams, data_psds(:,ii));
%     scatter(dsim, qsim(:,ii));
%     xlim([0, 20])
%     hold off
% end

%%  End the diary now because the program is finished

diary off
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

function prettyTimerDisplay(timeInSeconds)
    %   Remove whole number of hours first
    h = floor(timeInSeconds/3600);
    partialSeconds = timeInSeconds - h*3600;
    %   Remove whole number of minutes
    m = floor(partialSeconds/60);
    partialSeconds = partialSeconds - m*60;
    %   Print
    fprintf("%d hours, %d minutes, %d seconds", h,m,round(partialSeconds));
end