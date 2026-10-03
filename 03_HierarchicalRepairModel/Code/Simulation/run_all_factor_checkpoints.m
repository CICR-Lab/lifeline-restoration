function plan=run_all_factor_checkpoints(cfg,design)

plan=all_factor_checkpoint_plan(cfg,design);
if ~exist(cfg.paths.checkpointDir,'dir'), mkdir(cfg.paths.checkpointDir); end

useParallel=cfg.useParallel && license('test','Distrib_Computing_Toolbox');
if useParallel && isempty(gcp('nocreate'))
    parpool(cfg.maxWorkers);
end

schema=all_factor_result_schema();
for p=1:height(plan)
    spec=plan(p,:);
    finalPath=fullfile(cfg.paths.checkpointDir,spec.fileName{1});
    if valid_checkpoint(finalPath,cfg,spec,schema)
        fprintf('All-factor checkpoint %d already valid; skipping.\n', ...
            spec.checkpointID);
        continue;
    end

    caseIDs=spec.caseIDs{1};
    caseRows=design.cases(caseIDs,:);
    nJobs=height(caseRows)*cfg.nSeeds;
    resultMatrix=nan(nJobs,schema.nColumns);
    errorMessages=cell(nJobs,1);
    fprintf(['Running checkpoint %d: type %d, scramble %d, bases ' ...
        '%d-%d (%d trajectories).\n'],spec.checkpointID, ...
        spec.designTypeCode,spec.scrambleID,spec.firstBaseID, ...
        spec.lastBaseID,spec.expectedRows);

    if useParallel
        parfor job=1:nJobs
            [caseLocal,seedID]=decode_job(job,cfg.nSeeds);
            [resultMatrix(job,:),errorMessages{job}]=run_all_factor_case( ...
                caseRows(caseLocal,:),seedID,cfg,schema);
        end
    else
        for job=1:nJobs
            [caseLocal,seedID]=decode_job(job,cfg.nSeeds);
            [resultMatrix(job,:),errorMessages{job}]=run_all_factor_case( ...
                caseRows(caseLocal,:),seedID,cfg,schema);
        end
    end

    failed=resultMatrix(:,schema.index.status_code)~=0;
    if any(failed)
        reportPath=write_failure_report(cfg,spec,resultMatrix, ...
            errorMessages,schema,failed);
        firstFailed=find(failed,1,'first');
        error('RepairTask:AllFactorCheckpointFailed', ...
            ['Checkpoint %d produced %d failed trajectories. No ' ...
             'checkpoint was published. First failure: %s. See %s.'], ...
            spec.checkpointID,nnz(failed),errorMessages{firstFailed},reportPath);
    end

    resultNames=schema.names;
    checkpointMeta=make_meta(cfg,spec);
    save_checkpoint(finalPath,resultMatrix,resultNames,checkpointMeta, ...
        cfg,spec,schema);
end
end

function [caseLocal,seedID]=decode_job(job,nSeeds)
seedID=mod(job-1,nSeeds)+1;
caseLocal=floor((job-1)/nSeeds)+1;
end

function meta=make_meta(cfg,spec)
meta.modelVersion=cfg.modelVersion;
meta.experimentID=cfg.experimentID;
meta.randomGenerator=cfg.randomGenerator;
meta.checkpointID=spec.checkpointID;
meta.expectedRows=spec.expectedRows;
meta.caseIDs=spec.caseIDs{1};
end

function save_checkpoint(finalPath,resultMatrix,resultNames,checkpointMeta, ...
    cfg,spec,schema)
partialPath=[finalPath,'.partial.mat'];
delete_if_exists(partialPath);
partialCleanup=onCleanup(@() delete_if_exists(partialPath));
save(partialPath,'resultMatrix','resultNames','checkpointMeta','-v7.3');
if ~valid_checkpoint(partialPath,cfg,spec,schema)
    error('RepairTask:AllFactorCheckpointValidation', ...
        'New checkpoint failed validation: %s',partialPath);
end

delete_if_exists(finalPath);
[ok,message]=movefile(partialPath,finalPath,'f');
if ~ok, error('RepairTask:AllFactorCheckpointMove','%s',message); end
clear partialCleanup;
end

function ok=valid_checkpoint(path,cfg,spec,schema)
ok=false;
if ~exist(path,'file'), return; end
try
    x=load(path,'resultMatrix','resultNames','checkpointMeta');
    m=x.resultMatrix;
    names=x.resultNames(:);
    meta=x.checkpointMeta;
    if ~checkpoint_matches_experiment(meta,cfg) || ...
            ~strcmp(meta.randomGenerator,cfg.randomGenerator) || ...
            meta.checkpointID~=spec.checkpointID || ...
            meta.expectedRows~=spec.expectedRows || ...
            ~isequal(meta.caseIDs(:),spec.caseIDs{1}(:)) || ...
            size(m,1)~=spec.expectedRows || ...
            size(m,2)~=schema.nColumns || ~isequal(names,schema.names(:))
        return;
    end

    runID=m(:,schema.index.run_id);
    caseID=m(:,schema.index.case_id);
    seedID=m(:,schema.index.seed_id);
    if any(~isfinite(runID)) || any(runID~=floor(runID)) || ...
            any(caseID~=floor(caseID)) || any(seedID~=floor(seedID)) || ...
            any(~ismember(caseID,spec.caseIDs{1})) || ...
            any(seedID<1 | seedID>cfg.nSeeds) || ...
            any(runID~=(caseID-1).*cfg.nSeeds+seedID)
        return;
    end
    [c,s]=ndgrid(spec.caseIDs{1},1:cfg.nSeeds);
    if ~isequal(sort(runID),sort((c(:)-1).*cfg.nSeeds+s(:))), return; end

    status=m(:,schema.index.status_code);
    ratios=[schema.index.tau80_tau50,schema.index.tau90_tau80, ...
        schema.index.tau95_tau90,schema.index.tau80_tau50_step, ...
        schema.index.tau90_tau80_step,schema.index.tau95_tau90_step];
    kappa=[schema.index.kappa,schema.index.kappa_step];
    eta=[schema.index.eta90,schema.index.eta90_step];
    entropy=m(:,schema.index.task_layer_entropy_normalized);
    if any(status~=0) || any(~isfinite(m(:,[ratios,kappa,eta])),'all') || ...
            any(m(:,ratios)<1-cfg.numericTolerance,'all') || ...
            any(m(:,kappa)<=0,'all') || ...
            any(m(:,eta)<0 | m(:,eta)>1,'all') || ...
            any(~isfinite(entropy)) || ...
            any(entropy< -cfg.numericTolerance | entropy>1+cfg.numericTolerance) || ...
            any(m(:,schema.index.is_monotone)~=1) || ...
            any(m(:,schema.index.all_reconnected)~=1) || ...
            any(abs(m(:,schema.index.final_service)-1)>100*cfg.numericTolerance)
        return;
    end

    ok=true;
catch
    ok=false;
end
end

function reportPath=write_failure_report( ...
    cfg,spec,resultMatrix,errorMessages,schema,failed)
if ~exist(cfg.paths.failureDir,'dir'), mkdir(cfg.paths.failureDir); end
run_id=resultMatrix(failed,schema.index.run_id);
case_id=resultMatrix(failed,schema.index.case_id);
seed_id=resultMatrix(failed,schema.index.seed_id);
status_code=resultMatrix(failed,schema.index.status_code);
errorMessages=errorMessages(:);
error_message=errorMessages(failed);
report=table(run_id,case_id,seed_id,status_code,error_message);
reportPath=fullfile(cfg.paths.failureDir, ...
    sprintf('fail_cp_%05d.csv',spec.checkpointID));
writetable(report,reportPath);
end

function delete_if_exists(path)
if exist(path,'file'), delete(path); end
end
