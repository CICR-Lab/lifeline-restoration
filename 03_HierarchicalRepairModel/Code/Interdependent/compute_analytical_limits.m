function analytical = compute_analytical_limits(cfg)

theta = cfg.analyticalTheta(:);
piVec = cfg.coupled.piB(:)';
L = numel(piVec)-1;
model=struct('meanDuration',1,'releaseSpec',struct('type','identity'), ...
    'pi',piVec);
response=analytical_response('frontier',theta,0:L,model);
Rscarce=response.systemRestoredFraction;
model.distName='exponential'; model.durationCV=1;
response=analytical_response('parallel',theta,0:L,model);
Rparallel=response.systemRestoredFraction;

analytical(1).id = 'analytical_scarce';
analytical(1).label = 'Analytical - scarce capacity';
analytical(1).time = theta;
analytical(1).R = Rscarce;
analytical(1).metrics = restoration_curve_metrics(theta,Rscarce,'linear');

analytical(2).id = 'analytical_parallel';
analytical(2).label = 'Analytical - fully parallel';
analytical(2).time = theta;
analytical(2).R = Rparallel;
analytical(2).metrics = restoration_curve_metrics(theta,Rparallel,'linear');
end
