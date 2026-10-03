function [row,errorMessage]=run_all_factor_case(caseRow,seedID,cfg,schema)

row=base_row(caseRow,seedID,cfg,schema);
errorMessage='';
try
    [net,realized]=af_generate_repair_task_network(caseRow,seedID,cfg);
catch ME
    row(schema.index.status_code)=1;
    errorMessage=sprintf('%s: %s',ME.identifier,ME.message);
    return;
end

nCrews=max(1,min(net.N,round(caseRow.crewFraction.*net.N)));
row(schema.index.n_tasks_realized)=realized.N;
row(schema.index.crew_fraction_realized)=nCrews./net.N;
row(schema.index.task_layer_entropy_normalized)= ...
    realized.taskLayerEntropyNormalized;
row(schema.index.duration_depth_ratio_realized)= ...
    realized.durationDepthRatio;
row(schema.index.parent_outdegree_cv_realized)= ...
    realized.parentOutdegreeCV;
row(schema.index.min_layer_count)=realized.minLayerCount;
row(schema.index.effective_service_task_count)= ...
    realized.effectiveServiceTaskCount;
row(schema.index.parent_use_effective_fraction)= ...
    realized.parentUseEffectiveFraction;
row(schema.index.cv_duration_within_layer_realized)= ...
    realized.cvDurationWithinLayer;
row(schema.index.cv_weight_within_layer_realized)= ...
    realized.cvWeightWithinLayer;
row(schema.index.rho_dw_realized)=realized.rhoDurationWeight;

try
    [trajectory,diagnostic]=simulate_repair_network_qf( ...
        net,nCrews,caseRow.qF,cfg);
    metrics=compute_restoration_indicators(trajectory,cfg);
    row(schema.index.chi_frontier)=diagnostic.chiFrontier;
    row(schema.index.crew_utilization)=diagnostic.crewUtilization;
    row(schema.index.max_restoration_jump)=diagnostic.maxRestorationJump;
    row=store_metrics(row,metrics,diagnostic,schema);
    row(schema.index.status_code)=0;
catch ME
    row(schema.index.status_code)=2;
    errorMessage=sprintf('%s: %s',ME.identifier,ME.message);
end
end

function row=base_row(caseRow,seedID,cfg,schema)
row=nan(1,schema.nColumns);
batchID=ceil(seedID./cfg.seedBatchSize);
row(schema.index.run_id)=(caseRow.caseID-1).*cfg.nSeeds+seedID;
row(schema.index.case_id)=caseRow.caseID;
row(schema.index.design_type_code)=caseRow.designTypeCode;
row(schema.index.scramble_id)=caseRow.scrambleID;
row(schema.index.base_id_within_scramble)=caseRow.baseIDWithinScramble;
row(schema.index.base_id_global)=caseRow.baseIDGlobal;
row(schema.index.background_id)=caseRow.backgroundID;
row(schema.index.random_block_id)=caseRow.randomBlockID;
row(schema.index.component_block_id)=af_make_random_block_id( ...
    caseRow.randomBlockID,seedID,cfg);
row(schema.index.matrix_id)=caseRow.matrixID;
row(schema.index.hybrid_factor_id)=caseRow.hybridFactorID;
row(schema.index.hybrid_group_id)=caseRow.hybridGroupID;
row(schema.index.response_cell_id)=caseRow.responseCellID;
row(schema.index.response_axis1_id)=caseRow.responseAxis1ID;
row(schema.index.response_axis2_id)=caseRow.responseAxis2ID;
row(schema.index.seed_id)=seedID;
row(schema.index.seed_batch_id)=batchID;
row(schema.index.seed_within_batch)= ...
    seedID-(batchID-1).*cfg.seedBatchSize;
row(schema.index.n_tasks_target)=caseRow.nTasks;
row(schema.index.maximum_depth)=caseRow.maximumDepth;
row(schema.index.task_count_concentration)=caseRow.taskCountConcentration;
row(schema.index.parent_concentration)=caseRow.parentConcentration;
row(schema.index.duration_depth_ratio)=caseRow.durationDepthRatio;
row(schema.index.cv_duration_within_layer_target)= ...
    caseRow.cvDurationWithinLayer;
row(schema.index.cv_weight_target)=caseRow.cvWeightWithinLayer;
row(schema.index.rho_dw_target)=caseRow.rhoDurationWeight;
row(schema.index.crew_fraction_target)=caseRow.crewFraction;
row(schema.index.q_f_target)=caseRow.qF;
row(schema.index.status_code)=9;
end

function row=store_metrics(row,m,d,schema)
row(schema.index.tau80_tau50)=m.tau80_tau50;
row(schema.index.tau90_tau80)=m.tau90_tau80;
row(schema.index.tau95_tau90)=m.tau95_tau90;
row(schema.index.kappa)=m.kappa;
row(schema.index.eta90)=m.eta90;
row(schema.index.tau80_tau50_step)=m.step.tau80_tau50;
row(schema.index.tau90_tau80_step)=m.step.tau90_tau80;
row(schema.index.tau95_tau90_step)=m.step.tau95_tau90;
row(schema.index.kappa_step)=m.step.kappa;
row(schema.index.eta90_step)=m.step.eta90;
row(schema.index.final_service)=d.finalService;
row(schema.index.is_monotone)=d.isMonotone;
row(schema.index.all_reconnected)=d.allReconnected;
end
