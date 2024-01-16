function dy = NanoSimRHS(~,y,settings) %#codegen
    %   Chemical reactions
    dy = sum( -settings.A .* ( prod( y .^ max(0,settings.A) ) .* settings.K ) ,2);

    %   Particle Growth
    pstart = settings.pstart;
    pend = settings.pend;
    dGrowth = y(settings.gidx) * y(pstart:pend-1) .* settings.gKernel';
    dy(pstart:pend-1) = dy(pstart:pend-1) - dGrowth;
    dy(pstart+1:pend) = dy(pstart+1:pend) + dGrowth;
    dy(settings.gRxnIdx) = dy(settings.gRxnIdx) - settings.gRxnCoeff' .* sum(dGrowth);

    %   Particle Agglomeration
    for iii=pstart:settings.cutoff
        sizeiii = iii-1; %fixme
        jjj = iii:settings.cutoff;
        sizejjj = jjj-1; %fixme
        created = sizejjj + sizeiii + 1; %fixme
        dAggl = y(iii) * y(jjj) .* settings.aKernel(sizejjj,sizeiii);

        dy(iii) = dy(iii) - sum(dAggl);
        dy(jjj) = dy(jjj) - dAggl;
        dy(created) = dy(created) + dAggl;
    end
end









