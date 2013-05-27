% Copyright (C) Daphne Koller, Stanford University, 2012

function [MEU OptimalDecisionRule] = OptimizeMEU( I )

  % Inputs: An influence diagram I with a single decision node and a single utility node.
  %         I.RandomFactors = list of factors for each random variable.  These are CPDs, with
  %              the child variable = D.var(1)
  %         I.DecisionFactors = factor for the decision node.
  %         I.UtilityFactors = list of factors representing conditional utilities.
  % Return value: the maximum expected utility of I and an optimal decision rule 
  % (represented again as a factor) that yields that expected utility.
  
  % We assume I has a single decision node.
  % You may assume that there is a unique optimal decision.
  D = I.DecisionFactors(1);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % YOUR CODE HERE...
  % 
  % Some other information that might be useful for some implementations
  % (note that there are multiple ways to implement this):
  % 1.  It is probably easiest to think of two cases - D has parents and D 
  %     has no parents.
  % 2.  You may find the Matlab/Octave function setdiff useful.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
  
  % D has no parents
  if length(D.var) == 1
      OptimalDecisionRule = D;
      MEU = sum(D.val .* I.UtilityFactors.val);
  % D has parents
  else
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

end
