function design = build_all_factor_design(cfg)

sobolCases=build_sobol_cases(cfg);
taskServiceCases=build_task_service_cases(cfg);
capacityQCases=build_capacity_q_cases(cfg);

cases=[sobolCases;taskServiceCases;capacityQCases];
cases=sortrows(cases,{'designTypeCode','scrambleID', ...
    'baseIDWithinScramble','matrixID','responseCellID'});
cases.caseID=(1:height(cases))';
cases=movevars(cases,'caseID','Before',1);

if height(cases)~=cfg.expectedCases
    error('RepairTask:AllFactorCaseCount', ...
        'Expected %d design cases but generated %d.', ...
        cfg.expectedCases,height(cases));
end
design.cases=cases;
design.sobolCases=cases(cases.designTypeCode==1,:);
design.taskServiceCases=cases(cases.designTypeCode==2,:);
design.capacityQCases=cases(cases.designTypeCode==3,:);
end

function cases=build_sobol_cases(cfg)
n=cfg.nBasePerScramble;
m=cfg.nSobolMatrices;
tables=cell(cfg.nScrambles,1);
for scrambleID=1:cfg.nScrambles
    points=af_scrambled_sobol_points(n,2*cfg.nFactors, ...
        cfg.masterSeed+1000+scrambleID,cfg.sobolSkip);
    A=points(:,1:cfg.nFactors);
    B=points(:,cfg.nFactors+(1:cfg.nFactors));
    unit=zeros(n*m,cfg.nFactors);
    matrixID=zeros(n*m,1);
    matrixRole=cell(n*m,1);
    hybridFactorID=zeros(n*m,1);
    hybridFactorName=repmat({''},n*m,1);
    hybridGroupID=zeros(n*m,1);
    hybridGroupName=repmat({''},n*m,1);
    baseIDWithinScramble=zeros(n*m,1);
    row=0;
    for b=1:n
        matrixPoints=cell(m,1);
        matrixPoints{1}=A(b,:);
        matrixPoints{2}=B(b,:);
        for factorID=1:cfg.nFactors
            value=A(b,:);
            value(factorID)=B(b,factorID);
            matrixPoints{2+factorID}=value;
        end
        for groupID=1:numel(cfg.groupNames)
            value=A(b,:);
            ids=cfg.groupFactorIDs{groupID};
            value(ids)=B(b,ids);
            matrixPoints{2+cfg.nFactors+groupID}=value;
        end
        for matrix=1:m
            row=row+1;
            unit(row,:)=matrixPoints{matrix};
            matrixID(row)=matrix;
            baseIDWithinScramble(row)=b;
            if matrix==1
                matrixRole{row}='A';
            elseif matrix==2
                matrixRole{row}='B';
            elseif matrix<=2+cfg.nFactors
                factorID=matrix-2;
                matrixRole{row}='AB_factor';
                hybridFactorID(row)=factorID;
                hybridFactorName{row}=cfg.factorNames{factorID};
            else
                groupID=matrix-(2+cfg.nFactors);
                matrixRole{row}='AB_group';
                hybridGroupID(row)=groupID;
                hybridGroupName{row}=cfg.groupNames{groupID};
            end
        end
    end
    nRows=n*m;
    designTypeCode=ones(nRows,1);
    designType=repmat({'sobol'},nRows,1);
    scrambleColumn=repmat(scrambleID,nRows,1);
    baseIDGlobal=(scrambleID-1)*n+baseIDWithinScramble;
    backgroundID=zeros(nRows,1);
    randomBlockID=baseIDGlobal;
    responseCellID=zeros(nRows,1);
    responseAxis1ID=zeros(nRows,1);
    responseAxis2ID=zeros(nRows,1);
    meta=table(designTypeCode,designType,scrambleColumn, ...
        baseIDWithinScramble,baseIDGlobal,backgroundID,randomBlockID, ...
        matrixID,matrixRole,hybridFactorID,hybridFactorName, ...
        hybridGroupID,hybridGroupName,responseCellID, ...
        responseAxis1ID,responseAxis2ID,'VariableNames', ...
        {'designTypeCode','designType','scrambleID', ...
         'baseIDWithinScramble','baseIDGlobal','backgroundID', ...
         'randomBlockID','matrixID','matrixRole','hybridFactorID', ...
         'hybridFactorName','hybridGroupID','hybridGroupName', ...
         'responseCellID','responseAxis1ID','responseAxis2ID'});
    tables{scrambleID}=assemble_cases(meta,unit,cfg);
end
cases=vertcat(tables{:});
end

function cases=build_task_service_cases(cfg)
nPer=cfg.nBackgroundsPerResponseScramble;
n=cfg.nResponseBackgrounds;
nCell=cfg.nResponseCells;
unit=zeros(n*nCell,cfg.nFactors);
backgroundID=zeros(n*nCell,1);
scrambleID=zeros(n*nCell,1);
baseIDWithinScramble=zeros(n*nCell,1);
responseCellID=zeros(n*nCell,1);
responseAxis1ID=zeros(n*nCell,1);
responseAxis2ID=zeros(n*nCell,1);
row=0;
for s=1:cfg.nResponseScrambles
    background=af_scrambled_sobol_points(nPer,cfg.nFactors, ...
        cfg.masterSeed+2000+s,cfg.sobolSkip+4096);
    for b=1:nPer
        globalBackground=(s-1)*nPer+b;
        for taskID=1:numel(cfg.taskLayerProfiles)
            for serviceID=1:numel(cfg.serviceLayerProfiles)
                row=row+1;
                unit(row,:)=background(b,:);
                unit(row,3)=(taskID-0.5)/numel(cfg.taskLayerProfiles);
                unit(row,7)=(serviceID-0.5)/numel(cfg.serviceLayerProfiles);
                backgroundID(row)=globalBackground;
                scrambleID(row)=s;
                baseIDWithinScramble(row)=b;
                responseAxis1ID(row)=taskID;
                responseAxis2ID(row)=serviceID;
                responseCellID(row)=(taskID-1)* ...
                    numel(cfg.serviceLayerProfiles)+serviceID;
            end
        end
    end
end
cases=response_meta_and_map(2,'response_task_service',unit, ...
    backgroundID,scrambleID,baseIDWithinScramble,responseCellID, ...
    responseAxis1ID,responseAxis2ID, ...
    cfg.nBase,cfg);
end

function cases=build_capacity_q_cases(cfg)
nPer=cfg.nBackgroundsPerResponseScramble;
n=cfg.nResponseBackgrounds;
nCell=cfg.nResponseCells;
unit=zeros(n*nCell,cfg.nFactors);
backgroundID=zeros(n*nCell,1);
scrambleID=zeros(n*nCell,1);
baseIDWithinScramble=zeros(n*nCell,1);
responseCellID=zeros(n*nCell,1);
responseAxis1ID=zeros(n*nCell,1);
responseAxis2ID=zeros(n*nCell,1);
row=0;
for s=1:cfg.nResponseScrambles
    background=af_scrambled_sobol_points(nPer,cfg.nFactors, ...
        cfg.masterSeed+3000+s,cfg.sobolSkip+8192);
    for b=1:nPer
        globalBackground=(s-1)*nPer+b;
        for capacityID=1:numel(cfg.capacityAnchors)
            for qID=1:numel(cfg.qFAnchors)
                row=row+1;
                unit(row,:)=background(b,:);
                unit(row,12)=log(cfg.capacityAnchors(capacityID)/ ...
                    cfg.crewFractionBounds(1))/log(cfg.crewFractionBounds(2)/ ...
                    cfg.crewFractionBounds(1));
                unit(row,13)=cfg.qFAnchors(qID);
                backgroundID(row)=globalBackground;
                scrambleID(row)=s;
                baseIDWithinScramble(row)=b;
                responseAxis1ID(row)=capacityID;
                responseAxis2ID(row)=qID;
                responseCellID(row)=(capacityID-1)* ...
                    numel(cfg.qFAnchors)+qID;
            end
        end
    end
end
cases=response_meta_and_map(3,'response_capacity_q',unit, ...
    backgroundID,scrambleID,baseIDWithinScramble,responseCellID, ...
    responseAxis1ID,responseAxis2ID, ...
    cfg.nBase+cfg.nResponseBackgrounds,cfg);
for row=1:height(cases)
    cases.crewFraction(row)=cfg.capacityAnchors( ...
        cases.responseAxis1ID(row));
    cases.qF(row)=cfg.qFAnchors(cases.responseAxis2ID(row));
end
end

function cases=response_meta_and_map(code,name,unit,backgroundID, ...
    scrambleID,baseIDWithinScramble,responseCellID,responseAxis1ID, ...
    responseAxis2ID,blockOffset,cfg)
nRows=size(unit,1);
designTypeCode=repmat(code,nRows,1);
designType=repmat({name},nRows,1);
baseIDGlobal=blockOffset+backgroundID;
randomBlockID=baseIDGlobal;
matrixID=zeros(nRows,1);
matrixRole=repmat({name},nRows,1);
hybridFactorID=zeros(nRows,1);
hybridFactorName=repmat({''},nRows,1);
hybridGroupID=zeros(nRows,1);
hybridGroupName=repmat({''},nRows,1);
meta=table(designTypeCode,designType,scrambleID, ...
    baseIDWithinScramble,baseIDGlobal,backgroundID,randomBlockID, ...
    matrixID,matrixRole,hybridFactorID,hybridFactorName, ...
    hybridGroupID,hybridGroupName,responseCellID, ...
    responseAxis1ID,responseAxis2ID);
cases=assemble_cases(meta,unit,cfg);
end

function cases=assemble_cases(meta,unit,cfg)
mapped=map_all_factor_inputs(unit,cfg);
cases=[meta,mapped];
for factorID=1:cfg.nFactors
    cases.(sprintf('unitPoint%d',factorID))=unit(:,factorID);
end
end
