% function [nll, grad] = InstanceNegLogLikelihood(X, y, theta, modelParams)
% returns the negative log-likelihood and its gradient, given a CRF with parameters theta,
% on data (X, y).
%
% Inputs:
% X            Data.                           (numCharacters x numImageFeatures matrix)
%              X(:,1) is all ones, i.e., it encodes the intercept/bias term.
% y            Data labels.                    (numCharacters x 1 vector)
% theta        CRF weights/parameters.         (numParams x 1 vector)
%              These are shared among the various singleton / pairwise features.
% modelParams  Struct with three fields:
%   .numHiddenStates     in our case, set to 26 (26 possible characters)
%   .numObservedStates   in our case, set to 2  (each pixel is either on or off)
%   .lambda              the regularization parameter lambda
%
% Outputs:
% nll          Negative log-likelihood of the data.    (scalar)
% grad         Gradient of nll with respect to theta   (numParams x 1 vector)
%
% Copyright (C) Daphne Koller, Stanford Univerity, 2012

function [nll, grad] = InstanceNegLogLikelihood(X, y, theta, modelParams)

% featureSet is a struct with two fields:
%    .numParams - the number of parameters in the CRF (this is not numImageFeatures
%                 nor numFeatures, because of parameter sharing)
%    .features  - an array comprising the features in the CRF.
%
% Each feature is a binary indicator variable, represented by a struct
% with three fields:
%    .var          - a vector containing the variables in the scope of this feature
%    .assignment   - the assignment that this indicator variable corresponds to
%    .paramIdx     - the index in theta that this feature corresponds to
%
% For example, if we have:
%
%   feature = struct('var', [2 3], 'assignment', [5 6], 'paramIdx', 8);
%
% then feature is an indicator function over X_2 and X_3, which takes on a value of 1
% if X_2 = 5 and X_3 = 6 (which would be 'e' and 'f'), and 0 otherwise.
% Its contribution to the log-likelihood would be theta(8) if it's 1, and 0 otherwise.
%
% If you're interested in the implementation details of CRFs,
% feel free to read through GenerateAllFeatures.m and the functions it calls!
% For the purposes of this assignment, though, you don't
% have to understand how this code works. (It's complicated.)

featureSet = GenerateAllFeatures(X, modelParams);

% Use the featureSet to calculate nll and grad.
% This is the main part of the assignment, and it is very tricky - be careful!
% You might want to code up your own numerical gradient checker to make sure
% your answers are correct.
%
% Hint: you can use CliqueTreeCalibrate to calculate logZ effectively.
%       We have halfway-modified CliqueTreeCalibrate; complete our implementation
%       if you want to use it to compute logZ.

nll = 0;
grad = zeros(size(theta));
%%%
% Your code here:
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
        factor_1.val(idx) = factor_1.val(idx) + theta(feature.paramIdx);
    elseif feature.var == 2
        idx = AssignmentToIndex(feature.assignment, factor_2.card);
        factor_2.val(idx) = factor_2.val(idx) + theta(feature.paramIdx);
    elseif feature.var == 3
        idx = AssignmentToIndex(feature.assignment, factor_3.card);
        factor_3.val(idx) = factor_3.val(idx) + theta(feature.paramIdx);
    elseif feature.var == [1 2]
        idx = AssignmentToIndex(feature.assignment, factor_12.card);
        factor_12.val(idx) = factor_12.val(idx) + theta(feature.paramIdx);
    elseif feature.var == [2 3]
        idx = AssignmentToIndex(feature.assignment, factor_23.card);
        factor_23.val(idx) = factor_23.val(idx) + theta(feature.paramIdx);
    else
        disp(sprintf('unknown feature'));
    end
end
factors = [factor_1, factor_2, factor_3, factor_12, factor_23];
for i=1:length(factors)
    factors(i).val = exp(factors(i).val);
end
%%
factors = [];
for i=1:N
    feature = featureSet.features(i);
    factor = struct ('var', [], 'card', [], 'val', []);
    factor.var = feature.var;
    factor.card = 26 * ones(1, length(factor.var));
    factor.val = zeros(1, prod(factor.card));
    idx = AssignmentToIndex(feature.assignment, factor.card);
    factor.val(idx) = theta(feature.paramIdx);
    factors = [factors, factor];
end
%%
P = CreateCliqueTree(factors,[]);
[P logZ] = CliqueTreeCalibrate(P, false);
f1 = FactorMarginalization(P.cliqueList(1), 2);
f2 = FactorMarginalization(P.cliqueList(1), 1);
f3 = FactorMarginalization(P.cliqueList(2), 2);
f12 = P.cliqueList(1);
f23 = P.cliqueList(2);
factors_2 = [f1, f2, f3, f12, f23];
for i=1:length(factors_2)    
    factors_2(i) = NormalizeFactorValues(factors_2(i));
end
t_featurecounts = zeros(1, length(theta));
t_modelfeaturecounts = zeros(1, length(theta));
for i=1:length(t_featurecounts)
    for j=1:N
        feature = featureSet.features(j);
        if feature.paramIdx == i
            if feature.assignment == y(feature.var)
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
nll = logZ + modelParams.lambda / 2 * sum(theta.^2) - sum(theta .* t_featurecounts);
grad = t_modelfeaturecounts - t_featurecounts + modelParams.lambda * theta;
end
