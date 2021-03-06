clear;
load 'Part2Sample';
featureSet = GenerateAllFeatures(sampleX, sampleModelParams);

nll = 0;
grad = zeros(size(sampleTheta));

N = length(featureSet.features(:));
factor_1 = struct ('var', [], 'card', [], 'val', []);
factor_1.var = 1; factor_1.card = 26; factor_1.val = zeros(1, prod(factor_1.card));
factor_2 = struct ('var', [], 'card', [], 'val', []);
factor_2.var = 2; factor_2.card = 26; factor_2.val = zeros(1, prod(factor_2.card));
factor_3 = struct ('var', [], 'card', [], 'val', []);
factor_3.var = 3; factor_3.card = 26; factor_3.val = zeros(1, prod(factor_3.card));
factor_12 = struct ('var', [], 'card', [], 'val', []);
factor_12.var = [1 2]; factor_12.card = [26 26]; factor_12.val = zeros(1, prod(factor_12.card));
facotr_23 = struct ('var', [], 'card', [], 'val', []);
factor_23.var = [2 3]; factor_23.card = [26 26]; factor_23.val = zeros(1, prod(factor_23.card));

for i=1:N
    feature = featureSet.features(i);
    if feature.var == 1
        idx = AssignmentToIndex(feature.assignment, factor_1.card);
        factor_1.val(idx) = factor_1.val(idx) + sampleTheta(feature.paramIdx);
    elseif feature.var == 2
        idx = AssignmentToIndex(feature.assignment, factor_2.card);
        factor_2.val(idx) = factor_2.val(idx) + sampleTheta(feature.paramIdx);
    elseif feature.var == 3
        idx = AssignmentToIndex(feature.assignment, factor_3.card);
        factor_3.val(idx) = factor_3.val(idx) + sampleTheta(feature.paramIdx);
    elseif feature.var == [1 2]
        idx = AssignmentToIndex(feature.assignment, factor_12.card);
        factor_12.val(idx) = factor_12.val(idx) + sampleTheta(feature.paramIdx);
    elseif feature.var == [2 3]
        idx = AssignmentToIndex(feature.assignment, factor_23.card);
        factor_23.val(idx) = factor_23.val(idx) + sampleTheta(feature.paramIdx);
    else
        disp(sprintf('unknown feature'));
    end
end

factors = [factor_1, factor_2, factor_3, factor_12, factor_23];
for i=1:length(factors)
    factors(i).val = exp(factors(i).val);
    factors(i) = NormalizeFactorValues(factors(i));
end
P = CreateCliqueTree(factors,[]);
logZ = 0;
[P logZ] = CliqueTreeCalibrate(P, false);
f1 = FactorMarginalization(P.cliqueList(1), 2);
f2 = FactorMarginalization(P.cliqueList(1), 1);
f3 = FactorMarginalization(P.cliqueList(2), 2);
f12 = P.cliqueList(1);
f23 = P.cliqueList(2);
factors_2 = [f1, f2, f3, f12, f23];
for i=1:length(factors)    
    factors_2(i) = NormalizeFactorValues(factors_2(i));
end
t_featurecounts = zeros(1, length(sampleTheta));
t_modelfeaturecounts = zeros(1, length(sampleTheta));
for i=1:length(t_featurecounts)
    disp(sprintf('|DEBUG|param [%d]',i));
    for j=1:N
        feature = featureSet.features(j);
        if feature.paramIdx == i
            if feature.assignment == sampleY(feature.var)
                t_featurecounts(i) = t_featurecounts(i) + 1;
            end
            
            if feature.var == 1
                t_idx = AssignmentToIndex(feature.assignment, factors_2(1).card);
                t_modelfeaturecounts(i) = t_modelfeaturecounts(i) + factors_2(1).val(t_idx);
            elseif feature.var == 2
                t_idx = AssignmentToIndex(feature.assignment, factors_2(2).card);
                t_modelfeaturecounts(i) = t_modelfeaturecounts(i) + factors_2(2).val(t_idx);
            elseif feature.var == 3
                t_idx = AssignmentToIndex(feature.assignment, factors(3).card);
                t_modelfeaturecounts(i) = t_modelfeaturecounts(i) + factors_2(3).val(t_idx);
            elseif feature.var == [1 2]
                t_idx = AssignmentToIndex(feature.assignment, factors_2(4).card);
                t_modelfeaturecounts(i) = t_modelfeaturecounts(i) + factors_2(4).val(t_idx);
            elseif feature.var == [2 3]
                t_idx = AssignmentToIndex(feature.assignment, factors_2(5).card);
                t_modelfeaturecounts(i) = t_modelfeaturecounts(i) + factors_2(5).val(t_idx);
            else
            end
        end

    end
end

nll = logZ + sampleModelParams.lambda / 2 * sum(sampleTheta.^2) - sum(sampleTheta .* t_featurecounts);
grad = t_modelfeaturecounts - t_featurecounts + sampleModelParams.lambda * sampleTheta ;
