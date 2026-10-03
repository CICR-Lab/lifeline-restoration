function u = af_scrambled_sobol_points(n,dimension,seed,skip)

if exist('sobolset','file')==0
    error('RepairTask:SobolGeneratorUnavailable', ...
        ['sobolset and scramble are required for the direct all-factor ' ...
         'design. Install/enable Statistics and Machine Learning Toolbox.']);
end
previous=rng;
cleanup=onCleanup(@() rng(previous));
rng(double(seed),'twister');
try
    pointSet=sobolset(dimension,'Skip',double(skip));
    pointSet=scramble(pointSet,'MatousekAffineOwen');
    u=net(pointSet,n);
catch ME
    error('RepairTask:SobolGeneratorUnavailable', ...
        'Could not construct the required scrambled Sobol design: %s', ...
        ME.message);
end
u=min(max(double(u),eps('double')),1-eps('double'));
clear cleanup;
end
