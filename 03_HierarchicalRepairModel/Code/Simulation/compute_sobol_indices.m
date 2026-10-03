function estimate = compute_sobol_indices(data)

A=double(data.A);
B=double(data.B);
H=double(data.hybrid);
validate_sobol_data(A,B,H);
nPair = size(A,1);
nTarget = size(H,2);
nOutcome = size(A,3);

firstNumerator = nan(nTarget,nOutcome);
totalNumerator = nan(nTarget,nOutcome);
firstOrder = nan(nTarget,nOutcome);
totalOrder = nan(nTarget,nOutcome);
parameterVariance = nan(1,nOutcome);
grandMean = nan(1,nOutcome);
varianceMethod=repmat({'cross_batch_noise_corrected'},1,nOutcome);
allowFallback=logical(data.allowNonpositiveVarianceFallback);

for o = 1:nOutcome
    a1 = A(:,1,o); a2 = A(:,2,o);
    b1 = B(:,1,o); b2 = B(:,2,o);
    x1 = [a1; b1];
    x2 = [a2; b2];
    valid = isfinite(x1) & isfinite(x2);
    if nnz(valid) < 4
        error('AllFactorGSA:InsufficientSobolPairs', ...
            'Outcome %d has fewer than four valid cross-batch pairs.',o);
    end
    m1 = mean(x1(valid));
    m2 = mean(x2(valid));
    v = mean((x1(valid)-m1).*(x2(valid)-m2));
    useFallback=~isfinite(v) || v<=0;
    if useFallback
        if ~allowFallback
            error('AllFactorGSA:NonpositiveParameterVariance', ...
                ['Cross-batch parameter variance is nonpositive for outcome ' ...
                 '%d. Increase the base design or the number of seeds.'],o);
        end
        aMean=0.5*(a1+a2); bMean=0.5*(b1+b2);
        proxy=[aMean;bMean]; validProxy=isfinite(proxy);
        proxyMean=mean(proxy(validProxy));
        v=mean((proxy(validProxy)-proxyMean).^2);
        if ~isfinite(v) || v<=0
            error('AllFactorGSA:NonpositiveFallbackVariance', ...
                'Uncorrected quick-run fallback variance is nonpositive.');
        end
        varianceMethod{o}='uncorrected_batch_average_fallback';
        m1=proxyMean; m2=proxyMean;
    end
    parameterVariance(o) = v;
    grandMean(o) = 0.5*(m1+m2);

    for k = 1:nTarget
        h1 = H(:,k,1,o);
        h2 = H(:,k,2,o);
        validFirst = isfinite(b1) & isfinite(b2) & ...
            isfinite(h1) & isfinite(h2);
        validTotal = isfinite(a1) & isfinite(a2) & ...
            isfinite(h1) & isfinite(h2);
        if nnz(validFirst) < 2 || nnz(validTotal) < 2
            error('AllFactorGSA:InsufficientHybridPairs', ...
                'Target %d, outcome %d has insufficient valid hybrids.',k,o);
        end

        if useFallback
            bMean=0.5*(b1(validFirst)+b2(validFirst));
            hMean=0.5*(h1(validFirst)+h2(validFirst));
            firstNumerator(k,o)=mean((bMean-m1).*(hMean-m1));
            aMean=0.5*(a1(validTotal)+a2(validTotal));
            hMean=0.5*(h1(validTotal)+h2(validTotal));
            totalNumerator(k,o)=0.5*mean((aMean-hMean).^2);
        else
            nFirst12 = mean((b1(validFirst)-m1).* ...
                (h2(validFirst)-m2));
            nFirst21 = mean((b2(validFirst)-m2).* ...
                (h1(validFirst)-m1));
            firstNumerator(k,o) = 0.5*(nFirst12+nFirst21);

            d1 = a1(validTotal)-h1(validTotal);
            d2 = a2(validTotal)-h2(validTotal);
            totalNumerator(k,o) = 0.5*mean(d1.*d2);
        end
        firstOrder(k,o) = firstNumerator(k,o)./v;
        totalOrder(k,o) = totalNumerator(k,o)./v;
    end
end

estimate.nPairs = nPair;
estimate.nTargets = nTarget;
estimate.nOutcomes = nOutcome;
estimate.parameterVariance = parameterVariance;
estimate.firstNumerator = firstNumerator;
estimate.totalNumerator = totalNumerator;
estimate.firstOrderRaw = firstOrder;
estimate.totalOrderRaw = totalOrder;
estimate.firstOrderDisplay = min(max(firstOrder,0),1);
estimate.totalOrderDisplay = min(max(totalOrder,0),1);
estimate.interactionGapRaw = totalOrder-firstOrder;
estimate.grandMean = grandMean;
estimate.varianceMethod=varianceMethod;
estimate.usedUncorrectedFallback=strcmp( ...
    varianceMethod,'uncorrected_batch_average_fallback');
estimate.allowNonpositiveVarianceFallback=allowFallback;

sv=reshape(double(data.stochasticVariance),1,[]);
if numel(sv)~=nOutcome || any(~isfinite(sv) | sv<0)
    error('AllFactorGSA:InvalidStochasticVariance', ...
        'stochasticVariance must contain one finite nonnegative value per outcome.');
end
totalVariance=max(parameterVariance,0)+sv;
estimate.stochasticVariance=sv;
estimate.totalOutputVariance=totalVariance;
estimate.parameterShare=max(parameterVariance,0)./totalVariance;
estimate.stochasticShare=sv./totalVariance;

estimate.targetNames=reshape(data.targetNames,1,[]);
estimate.outcomeNames=reshape(data.outcomeNames,1,[]);
estimate.scrambleID=data.scrambleID;
estimate.baseID=data.baseID;
end

function validate_sobol_data(A,B,H)
if size(A,2)~=2 || size(B,2)~=2 || ~isequal(size(A),size(B))
    error('AllFactorGSA:InvalidBatchMeanShape', ...
        'A and B must have identical sizes and exactly two seed batches.');
end
if size(H,1)~=size(A,1) || size(H,2)<1 || ...
        size(H,3)~=2 || size(H,4)~=size(A,3)
    error('AllFactorGSA:InvalidHybridShape', ...
        ['hybrid must be nPair-by-nTarget-by-2-by-nOutcome and match ' ...
         'the A/B pair and outcome dimensions.']);
end
end
