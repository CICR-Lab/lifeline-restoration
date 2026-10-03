function sumCdf = iid_sum_cdf_fft(baseCdf,maxCount)
baseCdf = min(max(cummax(baseCdf(:)),0),1);
nGrid = numel(baseCdf);
if maxCount < 1 || maxCount ~= floor(maxCount)
    error('maxCount must be a positive integer.');
end

baseMass = [baseCdf(1); diff(baseCdf)];
baseMass = max(baseMass,0);
nFft = 2^nextpow2(maxCount*(nGrid-1)+1);
massTransform = fft(baseMass,nFft);
sumCdf = zeros(nGrid,maxCount);
for count = 1:maxCount
    mass = real(ifft(massTransform.^count));
    mass = max(mass(1:nGrid),0);
    sumCdf(:,count) = min(max(cummax(cumsum(mass)),0),1);
end
sumCdf(:,1) = baseCdf;
end
