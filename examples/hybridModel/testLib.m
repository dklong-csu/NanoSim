% clear variables
close all
% clc

%%  Test Iridium plus growth
kb = 1.37e5;
kf = 5e-7*kb;
k1 = 7.69e4;
k2 = 1.4e4;
k3 = 7.15e3;
k4 = 1.74e3;
M = 111;
S = 11.3;
T = 10.0;

%%
%   A ->[kf*S*S] As + L
%   As + L ->[kb] A
%   A + 2As ->[k1] P3 + L
%   A + Pi ->[G(k2,k3,M)*r(i)] Pi+1 + L
%   Particles -- 3:2500
nRxns = 3;
nSpecies = 3;
maxsize = 30000;
particles = 3:maxsize;
pstart = nSpecies+1;
pend = pstart + length(particles) - 1;
vecSize = pend; % + eMoM later

%   Reaction Coefficients
A = zeros(vecSize, nRxns);
%       A -> As + L
A(1:3,1) = [1;-1;-1];
%       As + L -> A
A(1:3,2) = [-1;1;1];
%       A + 2As -> P + L
A(1:4,3) = [1;2;-1;-1];
% A = sparse(A);
%   Reaction rates
K = [kf*S*S, kb, k1];

%   Particle Growth
%       Species that instigates growth
gidx = 1; % A + Pi -> Pi+1 + L
%       Species involved in the reaction
gRxnIdx = [1,3];
gRxnCoeff = [1,-1]; % >0 reactant, <0 product
%       growth kernel
%       Do not lose from the final particle so that mass is preserved
gKernel = growthFcn(particles(1:end-1), k2, k3, M);

%   Particle Agglomeration
aKernel = agglomFcn(M,k4);



%%
mySettings = NanoSimSettings(vecSize, ...
        A=A,...
        K=K,...
        pstart=pstart,...
        pend=pend,...
        gidx=gidx,...
        gKernel=gKernel,...
        gRxnIdx=gRxnIdx,...
        gRxnCoeff=gRxnCoeff,...
        cutoff=M,...
        aKernel=aKernel);

y0 = zeros(vecSize,1);
y0(1) = 0.0012;



%%  Solve ODE



userhsmex = true;
usejac = true;
usejacmex = true;



if userhsmex
    codegen -report NanoSimRHS.m -args {0.0, y0, mySettings}
    F = @(t,y) NanoSimRHS_mex(t,y,mySettings);
else
    F = @(t,y) NanoSimRHS(t,y,mySettings);
end


if usejac
    if usejacmex
        codegen -report NanoSimJac.m -args {0.0, y0, mySettings}
        odeOpts = odeset('Stats','off',...
            'RelTol',1e-8,...
            'AbsTol',1e-14,...
            'Jacobian',@(t,y) NanoSimJac_mex(t,y,mySettings),...
            'NonNegative',1:length(y0));
    else
        odeOpts = odeset('Stats','off',...
            'RelTol',1e-8,...
            'AbsTol',1e-14,...
            'Jacobian',@(t,y) NanoSimJac(t,y,mySettings),...
            'NonNegative',1:length(y0));
    end
else
    odeOpts = odeset('Stats','off',...
        'RelTol',1e-8,...
        'AbsTol',1e-14,...
        'NonNegative',1:length(y0));
end

tic
[t,y] = ode15s(F, [0 T], y0, odeOpts);
toc




%%  Plots
figure
%   Plot A, As, L
hold on
plot(t, y(:,1), 'DisplayName',"A")
plot(t, y(:,2), 'DisplayName',"As")
plot(t, y(:,3), 'DisplayName',"L")
xlabel("Time / h")
ylabel("Conc / mol/L")
legend
hold off

%   Plot PSD
figure;
ax = gca;
h = animatedline('LineStyle','-',...
    'Marker','none');
axis([0,maxsize,0,Inf])
atoms=3:maxsize;
n_plotpts = 10;
pt = round(linspace(1,length(t),n_plotpts));
for ttt=pt
    clearpoints(h)
    addpoints(h, atoms, y(ttt,4:end));
    title(ax,sprintf("PSD at %.3f hours",t(ttt)))
    drawnow
    pause(5/n_plotpts)
end


%%  Helper functions
function r = rFcn(size)
    r = 8/3 * (size .^ (2/3) );
end

function G = growthFcn(s,k1,k2,M)
    G = k1 * ones(size(s)) .* (s <= M) + k2 * ones(size(s)) .* (s > M);
    G = rFcn(s) .* G;
end

function A = agglomFcn(M,k)
    sizes = 1:M;
    A = rFcn(sizes)' * rFcn(sizes) * k;
end