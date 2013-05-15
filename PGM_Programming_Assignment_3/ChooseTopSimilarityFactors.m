function factors = ChooseTopSimilarityFactors (allFactors, F)
% This function chooses the similarity factors with the highest similarity
% out of all the possibilities.
%
% Input:
%   allFactors: An array of all the similarity factors.
%   F: The number of factors to select.
%
% Output:
%   factors: The F factors out of allFactors for which the similarity score
%     is highest.
%
% Hint: Recall that the similarity score for two images will be in every
%   factor table entry (for those two images' factor) where they are
%   assigned the same character value.
%
% Copyright (C) Daphne Koller, Stanford University, 2012

% If there are fewer than F factors total, just return all of them.
if (length(allFactors) <= F)
    factors = allFactors;
    return;
end

% Your code here:
% factors = allFactors; %%% REMOVE THIS LINE
K = allFactors(1).card(1);
N = length(allFactors);
fact = repmat(struct('var', [], 'card', [], 'val', []), N, 1);
for i=1:N
    fact(i).idx = i;
    fact(i).val = allFactors(i).val(AssignmentToIndex([1, 1], [K, K]));
end

[unused, order] = sort([fact(:).val], 'descend');
sorted_fact = fact(order);

factors = repmat(struct('var', [], 'card', [], 'val', []), F, 1);
for i=1:F
    factors(i) = allFactors(sorted_fact(i).idx);
end

end