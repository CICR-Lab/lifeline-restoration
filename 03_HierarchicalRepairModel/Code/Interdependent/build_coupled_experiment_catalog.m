function catalog=build_coupled_experiment_catalog(cfg)
base=base_scenario(cfg); catalog=repmat(base,0,1);
for gamma=cfg.coupled.gammaGrid
    for gB=cfg.coupled.gBGrid
        p=base; p.gammaTarget=gamma; p.gBTarget=gB;
        p.nRep=cfg.nRepCoupledCore;
        p.id=sprintf('core_gam%s_gb%s',tag(gamma),tag(gB));
        p.label=sprintf('Core: gamma=%.2g, g_B=%.2g',gamma,gB);
        p.group='core grid'; catalog(end+1,1)=p;
    end
end
ref=base; ref.gammaTarget=1; ref.gBTarget=.5;
ref.nRep=cfg.nRepCoupledSensitivity;
for q=[.5 .8 .95]
    p=ref; p.qE=q; p.id=['sens_qE_' tag(q)];
    p.label=sprintf('Support threshold q_E=%.2g',q); p.group='support threshold';
    catalog(end+1,1)=p;
end
for location={'shallow','intermediate','deep'}
    p=ref; p.gateLocation=location{1}; p.id=['sens_gate_' location{1}];
    p.label=['Gate location: ' location{1}]; p.group='gate location';
    catalog(end+1,1)=p;
end
for crewFraction=[.01 1]
    p=ref; p.crewFractionE=crewFraction; p.crewFractionB=crewFraction;
    p.id=['sens_capacity_' tag(crewFraction)];
    p.label=sprintf('Relative repair capacity C/N = %.2g',crewFraction);
    p.group='repair capacity';
    catalog(end+1,1)=p;
end
p=ref; p.supportRepresentation='task_specific'; p.id='sens_task_specific';
p.label='Task-specific electricity-support links'; p.group='support representation';
catalog(end+1,1)=p;
p=ref; p.CVdB=2; p.id='sens_high_CVd';
p.label='Dependent duration heterogeneity CV_d = 2'; p.group='task heterogeneity';
catalog(end+1,1)=p;
p=ref; p.CVwB=2; p.id='sens_high_CVw';
p.label='Dependent service-weight heterogeneity CV_w = 2'; p.group='task heterogeneity';
catalog(end+1,1)=p;
p=ref; p.qFB=0; p.id='sens_qF_0';
p.label='Dependent frontier-priority probability q_F = 0'; p.group='repair organization';
catalog(end+1,1)=p;
p=ref; p.qFB=.5; p.id='sens_qF_0p5';
p.label='Dependent frontier-priority probability q_F = 0.5'; p.group='repair organization';
catalog(end+1,1)=p;
p=ref; p.NE=250; p.NB=250;
p.id='sens_small_network'; p.label='Smaller paired networks'; p.group='network size';
catalog(end+1,1)=p;
p=ref; p.dependencyModeB='two_alternative_parents';
p.id='sens_alternative_paths'; p.label='Alternative dependent reconnection paths';
p.group='dependency structure'; catalog(end+1,1)=p;
end

function s=base_scenario(cfg)
c=cfg.coupled; g=cfg.generator;
s=struct('id','','label','','group','', ...
    'NE',c.NE,'NB',c.NB, ...
    'crewFractionE',c.crewFractionE,'crewFractionB',c.crewFractionB, ...
    'CVdE',c.CVdE,'CVdB',c.CVdB,'CVwE',c.CVwE,'CVwB',c.CVwB, ...
    'rhoDWE',c.rhoDWE,'rhoDWB',c.rhoDWB, ...
    'qFE',c.qFE,'qFB',c.qFB, ...
    'dependencyModeE',c.dependencyModeE,'dependencyModeB',c.dependencyModeB, ...
    'qE',c.qE,'gateLocation',c.gateLocation, ...
    'supportRepresentation',c.supportRepresentation, ...
    'gammaTarget',1,'gBTarget',.5,'nRep',cfg.nRepCoupledCore, ...
    'maximumDepthE',g.maximumDepthE,'maximumDepthB',g.maximumDepthB, ...
    'taskLayerProfileE',g.taskLayerProfileE,'taskLayerProfileB',g.taskLayerProfileB, ...
    'taskCountConcentrationE',g.taskCountConcentrationE, ...
    'taskCountConcentrationB',g.taskCountConcentrationB, ...
    'parentConcentrationE',g.parentConcentrationE, ...
    'parentConcentrationB',g.parentConcentrationB, ...
    'serviceLayerProfileE',g.serviceLayerProfileE, ...
    'serviceLayerProfileB',g.serviceLayerProfileB, ...
    'durationDepthRatioE',g.durationDepthRatioE, ...
    'durationDepthRatioB',g.durationDepthRatioB);
end

function value=tag(x)
value=strrep(sprintf('%.3g',x),'.','p'); value=strrep(value,'-','m');
end
