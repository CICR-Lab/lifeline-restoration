function profiles = service_profiles(L)
if L ~= 3
    error('The named Figure 4d profiles are defined for L = 3.');
end
l = 0:L;

shallow = [0.65 0.25 0.08 0.02];
uniform = ones(1,L+1)/(L+1);
intermediate = [0.1 0.4 0.4 0.1];
tail=shallow(end:-1:1);
names = {'Shallow-weighted','Approximately uniform', ...
    'Intermediate-weighted','Deep-tail'};
displayNames = {'Shallow','Uniform','Intermediate','Deep-tail'};
P = [shallow; uniform; intermediate; tail];
for k = 1:size(P,1)
    profiles(k).name = names{k};
    profiles(k).displayName = displayNames{k};
    profiles(k).pi = P(k,:);
    profiles(k).effectiveDepth = sum(l .* P(k,:));
end
end
