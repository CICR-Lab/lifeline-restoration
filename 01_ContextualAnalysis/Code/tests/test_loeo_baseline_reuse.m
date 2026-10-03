function tests=test_loeo_baseline_reuse
tests=functiontests(localfunctions);
end
function [D,y,c,p,cfg,B]=fixture
D=table([11;12;13],[1;1;2],'VariableNames',{'record_id','eid'});
y=[0;1;2];c=[.1;.9;1.8];p=[.05;.95;1.9];cfg=struct('random_seed',17,'bootstrap_replicates',2000);
S=table(.02,.12,.2,.4,'VariableNames',{'delta_ci95_low','delta_ci95_high','r2_ci95_low','r2_ci95_high'});
B=struct('record_id',D.record_id,'eid',D.eid,'observed',y, ...
    'context_prediction',c,'plus_prediction',p,'summary',S,'config',cfg);
B.d0=struct('record_id',D.record_id,'eid',D.eid,'observed',y,'prediction',c,'summary',S);
end
function testPairedReuse(t)
[D,y,c,p,cfg,B]=fixture;[lo,hi]=ctx_reuse_loeo_baseline(D,y,c,p,cfg,B,"paired");
verifyEqual(t,[lo hi],[.02 .12]);
end
function testD0Reuse(t)
[D,y,c,~,cfg,B]=fixture;[lo,hi]=ctx_reuse_loeo_baseline(D,y,c,[],cfg,B,"d0");
verifyEqual(t,[lo hi],[.2 .4]);
end
function testRecordMismatch(t)
[D,y,c,p,cfg,B]=fixture;D.record_id(1)=99;
verifyError(t,@()ctx_reuse_loeo_baseline(D,y,c,p,cfg,B,"paired"),'ContextualLOEO:BaselineMismatch');
end
function testEarthquakeMismatch(t)
[D,y,c,p,cfg,B]=fixture;D.eid(1)=3;
verifyError(t,@()ctx_reuse_loeo_baseline(D,y,c,p,cfg,B,"paired"),'ContextualLOEO:BaselineMismatch');
end
function testResponseMismatch(t)
[D,y,c,p,cfg,B]=fixture;y(1)=.1;
verifyError(t,@()ctx_reuse_loeo_baseline(D,y,c,p,cfg,B,"paired"),'ContextualLOEO:BaselineMismatch');
end
function testPredictionMismatch(t)
[D,y,c,p,cfg,B]=fixture;p(1)=.2;
verifyError(t,@()ctx_reuse_loeo_baseline(D,y,c,p,cfg,B,"paired"),'ContextualLOEO:BaselineMismatch');
end
function testConfigurationMismatch(t)
[D,y,c,p,cfg,B]=fixture;cfg.random_seed=18;
verifyError(t,@()ctx_reuse_loeo_baseline(D,y,c,p,cfg,B,"paired"),'ContextualLOEO:BaselineMismatch');
end
function testNonfinitePrediction(t)
[D,y,c,p,cfg,B]=fixture;p(1)=nan;
verifyError(t,@()ctx_reuse_loeo_baseline(D,y,c,p,cfg,B,"paired"),'ContextualLOEO:BaselineMismatch');
end
