% Copyright (C) Daphne Koller, Stanford University, 2012

function [MEU OptimalDecisionRule] = OptimizeWithJointUtility( I )
% Inputs: An influence diagram I with a single decision node and one or more utility nodes.
%         I.RandomFactors = list of factors for each random variable.  These are CPDs, with
%              the child variable = D.var(1)
%         I.DecisionFactors = factor for the decision node.
%         I.UtilityFactors = list of factors representing conditional utilities.
% Return value: the maximum expected utility of I and an optimal decision rule
% (represented again as a factor) that yields that expected utility.
% You may assume that there is a unique optimal decision.

% This is similar to OptimizeMEU except that we must find a way to
% combine the multiple utility factors.  Note: This can be done with very
% little code.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% YOUR CODE HERE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IU = struct('var', [], 'card', [], 'val', []);
for i=1:length(I.UtilityFactors)
    IU = FactorSum(IU,I.UtilityFactors(i));
end

I.UtilityFactors = IU;
D = I.DecisionFactors(1);
EUF = CalculateExpectedUtilityFactor(I);
EUF = ReorderFactorVars(EUF, D.var);
decisionVar = D.var(1);
conditionVars = setdiff(EUF.var, decisionVar);
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

end
