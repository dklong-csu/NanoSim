function settings = NanoSimSettings(vecSize, options)
%NanoSimSettings settings to control particle simulation
    arguments
        vecSize {mustBeInteger} = 1
        options.A {mustBeReal} = zeros(vecSize,1)
        options.K {mustBeReal} = 0
        
        options.pstart {mustBeInteger,mustBePositive} = 1
        options.pend {mustBeInteger,mustBePositive} = 1
        options.gidx {mustBeScalarOrEmpty,mustBeInteger,mustBePositive} = []
        options.gKernel {mustBeReal} = []
        options.gRxnIdx {mustBeInteger} = [];
        options.gRxnCoeff {mustBeInteger} = [];

        options.cutoff {mustBeScalarOrEmpty, mustBeInteger} = [];
        options.aKernel {mustBeReal} = [];
    end

    settings.vecSize = vecSize;
    settings.A = options.A;
    settings.K = options.K;

    settings.pstart = options.pstart;
    settings.pend = options.pend;
    settings.gidx = options.gidx;
    settings.gKernel = options.gKernel;
    settings.gRxnIdx = options.gRxnIdx;
    settings.gRxnCoeff = options.gRxnCoeff;

    settings.cutoff = options.cutoff;
    settings.aKernel = options.aKernel;
end







% settings.A
% settings.K
% settings.pstart
% settings.pend
% settings.gidx
% settings.G
% settings.gRxnIdx
% settings.gRxnCoeff