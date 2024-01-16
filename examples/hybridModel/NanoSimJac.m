function dfdy = NanoSimJac(~,y,settings) %#codegen
    dfdy = zeros(settings.vecSize);

    %   Chemical reactions
    for iii=1:settings.vecSize
        noniii = ~ismember(1:numel(y),iii);
        %   Chemical reactions
        dfdy(:,iii) = sum(...
            -settings.A .* ...
            ( ...
            y(iii) .^ max(settings.A(iii,:)-1,0) ...
            .* prod( y(noniii) .^ max(0, settings.A(noniii,:)) ) ...
            .* settings.K ...
            .* max(0, settings.A(iii,:)) ...
            ), 2);

        %   Particle Growth
        %   A + Pi -> Pi+1 + L
        %       For A and Piii
        if iii==settings.gidx
            %   A * sum(Pi*G) -> sum(Pi*G)
            dGrowthdy = y(settings.pstart:settings.pend-1) .* settings.gKernel';
            dfdy(settings.pstart:settings.pend-1, iii) = dfdy(settings.pstart:settings.pend-1, iii) - dGrowthdy;
            dfdy(settings.pstart+1:settings.pend, iii) = dfdy(settings.pstart+1:settings.pend, iii) + dGrowthdy;
            dfdy(settings.gRxnIdx, iii) = dfdy(settings.gRxnIdx, iii) - settings.gRxnCoeff' .* sum(dGrowthdy);
        elseif (iii >= settings.pstart) && (iii < settings.pend)
            %   A * Piii*Giii -> A * Giii
            %   A (neg)
            %   Piii (neg)
            %   Piii+1 (pos)
            %   L (pos)
            dGrowthdy = y(settings.gidx) * settings.gKernel(iii-settings.pstart+1);
            dfdy(iii,iii) = dfdy(iii,iii) - dGrowthdy;
            dfdy(iii+1,iii) = dfdy(iii+1,iii) + dGrowthdy;
            dfdy(settings.gRxnIdx, iii) = dfdy(settings.gRxnIdx, iii) - settings.gRxnCoeff' * dGrowthdy;
        end

        %   Particle Agglomeration
        %   Piii + P* -> Piii+*
        %   Deriv wrt yiii
        %       Aggl(iii,*) * Piii * P* -> Aggl(iii,*) * P*
        %   dfiiidiii = -sum(...)
        %   dfxdiii = - Aggl(iii,x) * 
        if (iii >= settings.pstart) && (iii < settings.cutoff)
            jjj = [settings.pstart:iii-1, iii+1:settings.cutoff];   % don't include iii because equation is special
            sizeiii = iii-1; %fixme
            sizejjj = jjj-1; %fixme
            created = sizejjj + sizeiii + 1; %fixme
            dAggldy = y(jjj) .* settings.aKernel(sizejjj, sizeiii);
            
            dfdy(jjj,iii) = dfdy(jjj,iii) - dAggldy;
            dfdy(created,iii) = dfdy(created,iii) + dAggldy;
            dfdy(iii,iii) = dfdy(iii,iii) - sum(dAggldy);
    
            %   Piii + Piii -> P2iii
            %   Special because it's a square and thus has differnt deriv
            %   fiii -= 2 * Piii^2 * Agg(iii)
            %   f2iii += Piii^2 * Agg(iii)
            %   dfiiidiii -= 4 * Piii * Agg(iii)
            %   df2iiidiii += 2 * Piii * Agg(iii)
            dAggl2dy = 2 * y(iii) * settings.aKernel(sizeiii,sizeiii);
            dfdy(iii,iii) = dfdy(iii,iii) - 2 * dAggl2dy;
            iii2idx = sizeiii*2 + 1; %fixme
            dfdy(iii2idx, iii) = dfdy(iii2idx, iii) + dAggl2dy;
        end
    end

    dfdy = sparse(dfdy);
end