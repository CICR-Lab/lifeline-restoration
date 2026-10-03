function v = indicator_vector(metrics)
v = [metrics.tau80_tau50 metrics.tau90_tau80 ...
    metrics.tau95_tau90 metrics.kappa metrics.eta90];
end
