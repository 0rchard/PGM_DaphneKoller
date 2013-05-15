%COMPUTEEXACTMARGINALSBP Runs exact inference and returns the marginals
%over all the variables (if isMax == 0) or the max-marginals (if isMax == 1). 
%
%   M = COMPUTEEXACTMARGINALSBP(F, E, isMax) takes a list of factors F,
%   evidence E, and a flag isMax, runs exact inference and returns the
%   final marginals for the variables in the network. If isMax is 1, then
%   it runs exact MAP inference, otherwise exact inference (sum-prod).
%   It returns an array of size equal to the number of variables in the 
%   network where M(i) represents the ith variable and M(i).val represents 
%   the marginals of the ith variable. 
%
% Copyright (C) Daphne Koller, Stanford University, 2012


function M = ComputeExactMarginalsBP(F, E, isMax)

% initialization
% you should set it to the correct value in your code
M = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% YOUR CODE HERE
%
% Implement Exact and MAP Inference.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
debug = 0;
P = CreateCliqueTree(F, E);
P = CliqueTreeCalibrate(P, isMax);
P = P.cliqueList;
N = length(P);
var_list = unique([F(:).var]);
V = length(var_list);

M = repmat(struct('var', [], 'card', [], 'val', []), V, 1);
for n=1:V
    v_candidate = var_list(n);
    for m = 1:N
        factor = P(m);
        if ~isempty(find(factor.var == v_candidate))
            ve = setdiff(factor.var, v_candidate);
            if isMax == 1;
                factor = FactorMaxMarginalization(factor, ve);
                M(n) = factor;
            else
                factor = FactorMarginalization(factor, ve);
                M(n) = FactorNormalization(factor);
            end            
            break;
        end
    end
end

end
