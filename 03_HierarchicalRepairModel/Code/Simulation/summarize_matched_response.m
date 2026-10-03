function response = summarize_matched_response(cases,caseMean,outcomeInfo,cfg,module, ...
    nBootstrap,bootstrapSeed)

if nargin<6 || isempty(nBootstrap), nBootstrap=0; end
if nargin<7 || isempty(bootstrapSeed), bootstrapSeed=260924; end
if ~isscalar(nBootstrap) || nBootstrap<0 || nBootstrap~=floor(nBootstrap)
    error('AllFactorGSA:ResponseBootstrapCount', ...
        'nBootstrap must be a nonnegative integer.');
end

cases=sortrows(cases,{'scrambleID','baseIDWithinScramble','responseCellID'});
backgrounds=unique(cases.backgroundID,'stable');
cells=unique(cases.responseCellID,'stable');
nBackground=numel(backgrounds); nCell=numel(cells); nOutcome=height(outcomeInfo);
if height(cases)~=nBackground*nCell
    error('AllFactorGSA:ResponseDesignBalance', ...
        '%s is not a complete background-by-cell cross.',module);
end

responseValues=nan(nBackground,nCell,nOutcome);
caseMap=nan(nBackground,nCell);
backgroundScramble=zeros(nBackground,1);
backgroundWithinScramble=zeros(nBackground,1);
for b=1:nBackground
    backgroundRows=cases(cases.backgroundID==backgrounds(b),:);
    backgroundScramble(b)=backgroundRows.scrambleID(1);
    backgroundWithinScramble(b)=backgroundRows.baseIDWithinScramble(1);
    for c=1:nCell
        hit=cases.backgroundID==backgrounds(b) & cases.responseCellID==cells(c);
        if nnz(hit)~=1
            error('AllFactorGSA:ResponseDesignBalance', ...
                '%s has a missing or duplicate matched cell.',module);
        end
        caseID=cases.caseID(hit);
        caseMap(b,c)=caseID;
        responseValues(b,c,:)=reshape(caseMean(caseID,:),1,1,nOutcome);
    end
end
scrambles=unique(backgroundScramble,'stable');
if numel(scrambles)~=cfg.nResponseScrambles
    error('AllFactorGSA:ResponseScrambleCount', ...
        '%s has %d response scrambles; expected %d.',module, ...
        numel(scrambles),cfg.nResponseScrambles);
end
bootstrapMeans=response_scramble_bootstrap(responseValues, ...
    backgroundScramble,nBootstrap,bootstrapSeed);

nLong=nBackground*nCell*nOutcome;
moduleColumn=repmat({module},nLong,1);
backgroundID=zeros(nLong,1); scrambleID=zeros(nLong,1);
baseIDWithinScramble=zeros(nLong,1); responseCellID=zeros(nLong,1);
axis1ID=zeros(nLong,1); axis2ID=zeros(nLong,1);
axis1Value=nan(nLong,1); axis2Value=nan(nLong,1);
axis1Label=cell(nLong,1); axis2Label=cell(nLong,1);
outcomeName=cell(nLong,1); seedMean=nan(nLong,1); row=0;
for b=1:nBackground
    for c=1:nCell
        meta=cases(cases.caseID==caseMap(b,c),:);
        [v1,v2,l1,l2]=axis_metadata(meta,cfg,module);
        for o=1:nOutcome
            row=row+1; backgroundID(row)=backgrounds(b);
            scrambleID(row)=backgroundScramble(b);
            baseIDWithinScramble(row)=backgroundWithinScramble(b);
            responseCellID(row)=cells(c);
            axis1ID(row)=meta.responseAxis1ID; axis2ID(row)=meta.responseAxis2ID;
            axis1Value(row)=v1; axis2Value(row)=v2;
            axis1Label{row}=l1; axis2Label{row}=l2;
            outcomeName{row}=outcomeInfo.outcomeName{o};
            seedMean(row)=responseValues(b,c,o);
        end
    end
end
backgroundTable=table(moduleColumn,scrambleID,baseIDWithinScramble, ...
    backgroundID,responseCellID,axis1ID,axis2ID, ...
    axis1Value,axis2Value,axis1Label,axis2Label,outcomeName,seedMean, ...
    repmat(cfg.nSeeds,nLong,1),'VariableNames', ...
    {'module','scrambleID','baseIDWithinScramble','backgroundID', ...
    'responseCellID','responseAxis1ID', ...
    'responseAxis2ID','responseAxis1Value','responseAxis2Value', ...
    'responseAxis1Label','responseAxis2Label','outcomeName', ...
    'seedMean','nSeeds'});

nCellRows=nCell*nOutcome;
moduleColumn=repmat({module},nCellRows,1);
responseCellID=zeros(nCellRows,1); axis1ID=zeros(nCellRows,1); axis2ID=zeros(nCellRows,1);
axis1Value=nan(nCellRows,1); axis2Value=nan(nCellRows,1);
axis1Label=cell(nCellRows,1); axis2Label=cell(nCellRows,1);
outcomeName=cell(nCellRows,1); meanValue=nan(nCellRows,1);
medianValue=meanValue; p05=meanValue; p10=meanValue;
p90=meanValue; p95=meanValue; row=0;
meanP025=meanValue; meanP975=meanValue;
for c=1:nCell
    meta=cases(cases.responseCellID==cells(c),:); meta=meta(1,:);
    [v1,v2,l1,l2]=axis_metadata(meta,cfg,module);
    for o=1:nOutcome
        row=row+1; x=responseValues(:,c,o);
        responseCellID(row)=cells(c); axis1ID(row)=meta.responseAxis1ID;
        axis2ID(row)=meta.responseAxis2ID; axis1Value(row)=v1; axis2Value(row)=v2;
        axis1Label{row}=l1; axis2Label{row}=l2; outcomeName{row}=outcomeInfo.outcomeName{o};
        meanValue(row)=mean(x,'omitnan'); medianValue(row)=quantile_local(x,0.5);
        p05(row)=quantile_local(x,0.05); p10(row)=quantile_local(x,0.1);
        p90(row)=quantile_local(x,0.9); p95(row)=quantile_local(x,0.95);
        if nBootstrap>0
            meanP025(row)=quantile_local(bootstrapMeans(:,c,o),0.025);
            meanP975(row)=quantile_local(bootstrapMeans(:,c,o),0.975);
        end
    end
end
cellTable=table(moduleColumn,responseCellID,axis1ID,axis2ID,axis1Value,axis2Value, ...
    axis1Label,axis2Label,outcomeName,repmat(nBackground,nCellRows,1), ...
    repmat(cfg.nResponseScrambles,nCellRows,1), ...
    repmat(cfg.nBackgroundsPerResponseScramble,nCellRows,1), ...
    repmat(cfg.nSeeds,nCellRows,1),meanValue,meanP025,meanP975, ...
    medianValue,p05,p10,p90,p95,repmat(nBootstrap,nCellRows,1), ...
    repmat({'complete_response_scramble_block'},nCellRows,1), ...
    'VariableNames',{'module','responseCellID','responseAxis1ID', ...
    'responseAxis2ID','responseAxis1Value','responseAxis2Value', ...
    'responseAxis1Label','responseAxis2Label','outcomeName','nBackgrounds', ...
    'nResponseScrambles','nBackgroundsPerScramble','nSeeds', ...
    'mean','meanP025','meanP975','median','p05','p10','p90','p95', ...
    'nBootstrap','bootstrapUnit'});

nPair=nCell*(nCell-1)/2; nContrast=nPair*nOutcome;
moduleColumn=repmat({module},nContrast,1);
referenceCellID=zeros(nContrast,1); comparisonCellID=zeros(nContrast,1);
referenceAxis1ID=zeros(nContrast,1); referenceAxis2ID=zeros(nContrast,1);
comparisonAxis1ID=zeros(nContrast,1); comparisonAxis2ID=zeros(nContrast,1);
referenceAxis1Label=cell(nContrast,1); referenceAxis2Label=cell(nContrast,1);
comparisonAxis1Label=cell(nContrast,1); comparisonAxis2Label=cell(nContrast,1);
outcomeName=cell(nContrast,1); meanDifference=nan(nContrast,1);
medianDifference=meanDifference; p10Difference=meanDifference; p90Difference=meanDifference;
meanDifferenceP025=meanDifference; meanDifferenceP975=meanDifference;
row=0;
for r=1:nCell-1
    refMeta=cases(cases.responseCellID==cells(r),:); refMeta=refMeta(1,:);
    [~,~,rl1,rl2]=axis_metadata(refMeta,cfg,module);
    for c=r+1:nCell
        compMeta=cases(cases.responseCellID==cells(c),:); compMeta=compMeta(1,:);
        [~,~,cl1,cl2]=axis_metadata(compMeta,cfg,module);
        for o=1:nOutcome
            row=row+1; d=responseValues(:,c,o)-responseValues(:,r,o);
            referenceCellID(row)=cells(r); comparisonCellID(row)=cells(c);
            referenceAxis1ID(row)=refMeta.responseAxis1ID;
            referenceAxis2ID(row)=refMeta.responseAxis2ID;
            comparisonAxis1ID(row)=compMeta.responseAxis1ID;
            comparisonAxis2ID(row)=compMeta.responseAxis2ID;
            referenceAxis1Label{row}=rl1; referenceAxis2Label{row}=rl2;
            comparisonAxis1Label{row}=cl1; comparisonAxis2Label{row}=cl2;
            outcomeName{row}=outcomeInfo.outcomeName{o};
            meanDifference(row)=mean(d,'omitnan');
            medianDifference(row)=quantile_local(d,0.5);
            p10Difference(row)=quantile_local(d,0.1);
            p90Difference(row)=quantile_local(d,0.9);
            if nBootstrap>0
                bootDifference=bootstrapMeans(:,c,o)-bootstrapMeans(:,r,o);
                meanDifferenceP025(row)=quantile_local(bootDifference,0.025);
                meanDifferenceP975(row)=quantile_local(bootDifference,0.975);
            end
        end
    end
end
contrastTable=table(moduleColumn,referenceCellID,comparisonCellID, ...
    referenceAxis1ID,referenceAxis2ID,comparisonAxis1ID,comparisonAxis2ID, ...
    referenceAxis1Label,referenceAxis2Label,comparisonAxis1Label, ...
    comparisonAxis2Label,outcomeName,repmat(nBackground,nContrast,1), ...
    repmat(cfg.nResponseScrambles,nContrast,1), ...
    repmat(cfg.nBackgroundsPerResponseScramble,nContrast,1), ...
    repmat(cfg.nSeeds,nContrast,1),meanDifference,meanDifferenceP025, ...
    meanDifferenceP975,medianDifference,p10Difference,p90Difference, ...
    repmat(nBootstrap,nContrast,1), ...
    repmat({'complete_response_scramble_block'},nContrast,1), ...
    'VariableNames',{'module','referenceCellID','comparisonCellID', ...
    'referenceAxis1ID','referenceAxis2ID','comparisonAxis1ID', ...
    'comparisonAxis2ID','referenceAxis1Label','referenceAxis2Label', ...
    'comparisonAxis1Label','comparisonAxis2Label','outcomeName', ...
    'nBackgrounds','nResponseScrambles','nBackgroundsPerScramble', ...
    'nSeeds','meanDifference','meanDifferenceP025', ...
    'meanDifferenceP975','medianDifference','p10Difference', ...
    'p90Difference','nBootstrap','bootstrapUnit'});

response.backgroundTable=backgroundTable;
response.cellTable=cellTable;
response.contrastTable=contrastTable;
response.dispersionDefinition= ...
    'P05-P95 and P10-P90 across matched background seed means; not confidence intervals';
response.excludedFromSobol=true;
response.excludedFromGlobalConsistency=true;
response.bootstrap=struct('nBootstrap',nBootstrap, ...
    'bootstrapSeed',bootstrapSeed, ...
    'bootstrapUnit','complete_response_scramble_block');
end

function boot=response_scramble_bootstrap(values,blockID,nBootstrap,bootstrapSeed)
nBackground=size(values,1); nCell=size(values,2); nOutcome=size(values,3);
boot=nan(nBootstrap,nCell,nOutcome);
if nBootstrap==0, return; end
blocks=unique(blockID,'stable');
members=cell(numel(blocks),1);
for s=1:numel(blocks)
    members{s}=find(blockID==blocks(s));
end
blockSize=cellfun(@numel,members);
if any(blockSize~=blockSize(1)) || sum(blockSize)~=nBackground
    error('AllFactorGSA:ResponseScrambleBalance', ...
        'Response scramble blocks must have equal background counts.');
end
stream=RandStream('mt19937ar','Seed',bootstrapSeed);
for b=1:nBootstrap
    selectedBlocks=randi(stream,numel(blocks),[numel(blocks),1]);
    selected=vertcat(members{selectedBlocks});
    boot(b,:,:)=reshape(mean(values(selected,:,:),1,'omitnan'), ...
        1,nCell,nOutcome);
end
end

function [v1,v2,l1,l2]=axis_metadata(meta,cfg,module)
i=meta.responseAxis1ID; j=meta.responseAxis2ID;
switch module
    case 'task_service'
        v1=i; v2=j; l1=cfg.taskLayerProfiles{i}; l2=cfg.serviceLayerProfiles{j};
    case 'capacity_q'
        v1=cfg.capacityAnchors(i); v2=cfg.qFAnchors(j);
        l1=sprintf('%.6g',v1); l2=sprintf('q_F=%.6g',v2);
    otherwise
        error('AllFactorGSA:UnknownResponseModule','Unknown response module %s.',module);
end
end

