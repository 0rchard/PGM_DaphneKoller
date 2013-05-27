% Copyright (C) Daphne Koller, Stanford University, 2012

function [MEU OptimalDecisionRule] = OptimizeLinearExpectations( I )
% Inputs: An influence diagram I with a single decision node and one or more utility nodes.
%         I.RandomFactors = list of factors for each random variable.  These are CPDs, with
%              the child variable = D.var(1)
%         I.DecisionFactors = factor for the decision node.
%         I.UtilityFactors = list of factors representing conditional utilities.
% Return value: the maximum expected utility of I and an optimal decision rule
% (represented again as a factor) that yields that expected utility.
% You may assume that there is a unique optimal decision.
%
% This is similar to OptimizeMEU except that we will have to account for
% multiple utility factors.  We will do this by calculating the expected
% utility factors and combining them, then optimizing with respect to that
% combined expected utility factor.
MEU = [];
OptimalDecisionRule = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% YOUR CODE HERE
%
% A decision rule for D assigns, for each joint assignment to D's parents,
% probability 1 to the best option from the EUF for that joint assignment
% to D's parents, and 0 otherwise.  Note that when D has no parents, it is
% a degenerate case we can handle separately for convenience.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = length(I.UtilityFactors);
EUF = struct('var', [], 'card', [], 'val', []);

% todo
for i=1:N
    I_t = I;
    I_t.UtilityFactors = I.UtilityFactors(i);
    EUF_t = CalculateExpectedUtilityFactor(I_t);
    for j=1:length(EUF_t)
        EUF = FactorSum(EUF, EUF_t(j));
    end
end

D = I.DecisionFactors;
EUF = ReorderFactorVars(EUF, D.var);
decisionVar = D.var(1);

OptimalDecisionRule = D;
OptimalDecisionRule.val = zeros(1, prod(D.card));
for i=1:length(EUF.val)
    assign = IndexToAssignment(i, EUF.card);
    conditionVarAssign = assign(2:length(assign));
    val_group = zeros(1, EUF.card(1));
    for j=1:EUF.card(1)
        t_assign = [j conditionVarAssign];
        idx = AssignmentToIndex(t_assign, EUF.card);
        val_group(j) = EUF.val(idx);
    end
    [value, idx] = max(val_group);
    t_assign = [idx conditionVarAssign];
    OptimalDecisionRule.val(AssignmentToIndex(t_assign, EUF.card)) = 1;
end
MEU = sum(OptimalDecisionRule.val .* EUF.val);