classdef ContextualInferenceTest < matlab.unittest.TestCase
    methods (TestClassSetup)
        function addImplementation(testCase)
            codeDir = fileparts(fileparts(mfilename("fullpath")));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(codeDir));
        end
    end

    methods (Test)
        function existingInvarianceChecks(testCase)
            previous = rng;
            testCase.addTeardown(@() rng(previous));
            result = validate_ctx_cr2_satterthwaite();
            testCase.verifyEqual(string(result.status), "complete");
        end

        function balancedClustersMatchMeanTest(testCase)
            clusterMeans = [1; 2; 4; 5; 8; 10];
            y = repelem(clusterMeans, 3) + repmat([-0.2; 0; 0.2], 6, 1);
            result = ctx_cr2_satterthwaite(ones(18,1), y, repelem((1:6)',3), 1, "mean");
            se = std(clusterMeans)/sqrt(6);
            testCase.verifyEqual(result.estimate, mean(clusterMeans), AbsTol=1e-12);
            testCase.verifyEqual(result.se_cr2, se, AbsTol=1e-12);
            testCase.verifyEqual(result.satterthwaite_df, 5, AbsTol=1e-12);
            testCase.verifyEqual(result.ci95_low, mean(clusterMeans)-tinv(.975,5)*se, AbsTol=1e-12);
            testCase.verifyEqual(result.p_value, 2*tcdf(-abs(mean(clusterMeans)/se),5), AbsTol=1e-12);
        end

        function missingRowsFail(testCase)
            testCase.verifyError(@() ctx_cr2_satterthwaite(ones(5,1),(1:4)',(1:5)'), ...
                "ContextualCR2:RowMismatch");
        end

        function rankDeficiencyFails(testCase)
            testCase.verifyError(@() ctx_cr2_satterthwaite(ones(5,2),(1:5)',(1:5)'), ...
                "ContextualCR2:RankDeficient");
        end

        function singularClusterFails(testCase)
            testCase.verifyError(@() ctx_cr2_satterthwaite(ones(5,1),(1:5)',ones(5,1)), ...
                "ContextualCR2:SingularAdjustment");
        end
    end
end
