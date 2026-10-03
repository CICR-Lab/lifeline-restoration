function plan = all_factor_checkpoint_plan(cfg,design)

rows=cell(0,1);
checkpointID=0;
for designTypeCode=1:3
    subset=design.cases(design.cases.designTypeCode==designTypeCode,:);
    for scrambleID=unique(subset.scrambleID(:)')
        local=subset(subset.scrambleID==scrambleID,:);
        baseIDs=unique(local.baseIDWithinScramble(:)');
        for firstIndex=1:cfg.basePointsPerCheckpoint:numel(baseIDs)
            checkpointID=checkpointID+1;
            selected=baseIDs(firstIndex:min( ...
                firstIndex+cfg.basePointsPerCheckpoint-1, ...
                numel(baseIDs)));
            mask=ismember(local.baseIDWithinScramble,selected);
            caseIDs=local.caseID(mask);
            expectedCases=numel(caseIDs);
            expectedRows=expectedCases*cfg.nSeeds;
            fileName=sprintf('cp_t%d_s%02d_b%05d-%05d.mat', ...
                designTypeCode,scrambleID,min(selected),max(selected));
            rows{checkpointID,1}={checkpointID,designTypeCode, ...
                scrambleID,min(selected),max(selected),expectedCases, ...
                expectedRows,fileName,caseIDs};
        end
    end
end
plan=cell2table(vertcat(rows{:}),'VariableNames', ...
    {'checkpointID','designTypeCode','scrambleID','firstBaseID', ...
     'lastBaseID','expectedCases','expectedRows','fileName','caseIDs'});
end
