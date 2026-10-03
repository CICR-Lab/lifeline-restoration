function signature = scenario_signature(scenario,algorithmVersion)

ignored = {'id','label','group','tags','nRep','seedKey'};
definition = scenario;
for k = 1:numel(ignored)
    if isfield(definition,ignored{k})
        definition = rmfield(definition,ignored{k});
    end
end
payload = jsonencode(struct('algorithmVersion',algorithmVersion, ...
    'scenario',orderfields(definition)));
signature = compute_sha256(payload);
end
