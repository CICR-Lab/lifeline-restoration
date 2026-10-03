function coupling=compute_coupling_analytical(analytical,cfg)

gamma=logspace(log10(0.2),log10(cfg.coupled.gammaPlotMax),241);
g=linspace(0,0.8,161);
[GAMMA,GB]=meshgrid(gamma,g);

for a=1:numel(analytical)
    t=analytical(a).time; RB=analytical(a).R;
    if a==1, limitName='frontier'; else, limitName='parallel'; end
    model=struct('meanDuration',1,'releaseSpec',struct('type','identity'), ...
        'pi',cfg.coupled.piE(:)');
    if strcmp(limitName,'parallel')
        model.distName='exponential'; model.durationCV=1;
    end
    response=analytical_response(limitName,t,0:numel(model.pi)-1,model);
    RE=response.systemRestoredFraction;
    uE=milestone_time(t,RE,cfg.coupled.qE,'linear');
    psi=trapz(t,1-RB);
    t90E=milestone_time(t,RE,0.9,'linear');
    t90B=analytical(a).metrics.tau90;
    EUB=(GB*uE)./(GAMMA*psi+GB*uE+eps);
    coupling.limit(a).id=analytical(a).id;
    coupling.limit(a).label=analytical(a).label;
    coupling.limit(a).uE=uE;
    coupling.limit(a).psi=psi;
    coupling.limit(a).rho90=GAMMA*(t90B/t90E);
    coupling.limit(a).epsilonUB=EUB;
end
coupling.gamma=gamma;
coupling.gB=g;
coupling.GAMMA=GAMMA;
coupling.GB=GB;
coupling.levels=cfg.coupled.mainContourLevels;
end
