load 'TestI0' TestI0;
OldD = TestI0.DecisionFactors;
TD = struct('var', [], 'card', [], 'val', []);
TD.var = [9 11];
TD.card = [2 2];
TD.val = [0 1 1 0];


% test_1
TestI0.DecisionFactors = OldD;
MEU_I = OptimizeWithJointUtility(TestI0);
TestI0.DecisionFactors = TD;
MEU_TD = OptimizeWithJointUtility(TestI0);

if MEU_TD > MEU_I
    disp(sprintf('%f,%f,%f',MEU_TD, MEU_I, 100 * log(MEU_TD - MEU_I)));
else
    disp(sprintf('%f,%f',MEU_TD, MEU_I));
end


% test_2
TestI0.RandomFactors(10).val = [0.999, 0.001, 0.25, 0.75];
TestI0.DecisionFactors = OldD;
MEU_I = OptimizeWithJointUtility(TestI0);
TestI0.DecisionFactors = TD;
MEU_TD = OptimizeWithJointUtility(TestI0);
if MEU_TD > MEU_I
    disp(sprintf('%f,%f,%f',MEU_TD, MEU_I, 100 * log(MEU_TD - MEU_I)));
else
    disp(sprintf('%f,%f',MEU_TD, MEU_I));
end

% test_3
TestI0.RandomFactors(10).val = [0.999, 0.001, 0.001, 0.999];
TestI0.DecisionFactors = OldD;
MEU_I = OptimizeWithJointUtility(TestI0);
TestI0.DecisionFactors = TD;
MEU_TD = OptimizeWithJointUtility(TestI0);
if MEU_TD > MEU_I
    disp(sprintf('%f,%f,%f',MEU_TD, MEU_I, 100 * log(MEU_TD - MEU_I)));
else
    disp(sprintf('%f,%f',MEU_TD, MEU_I));
end