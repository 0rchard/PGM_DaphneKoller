function factors = ComputeTripletFactors (images, tripletList, K)
% This function computes the triplet factor values for one word.
%
% Input:
%   images: An array of structs containing the 'img' value for each
%     character in the word.
%   tripletList: An array of the character triplets we will consider (other
%     factor values should be 1). tripletList(i).chars gives character
%     assignment, and triplistList(i).factorVal gives the value for that
%     entry in the factor table.
%   K: The alphabet size (accessible in imageModel.K for the provided
%     imageModel).
%
% Hint: Every character triple in the word will use the same 'val' table.
%   Consider computing that array once and then resusing for each factor.
%
% Copyright (C) Daphne Koller, Stanford University, 2012


n = length(images);

% If the word has fewer than three characters, then return an empty list.
if (n < 3)
    factors = [];
    return
end

factors = repmat(struct('var', [], 'card', [], 'val', []), n - 2, 1);

% Your code here:
valList = ones(1, K*K*K);
card = [K,K,K];
for i=1:length(tripletList)
    chars = tripletList(i).chars;
    c1 = chars(1);
    c2 = chars(2);
    c3 = chars(3);
    idx = AssignmentToIndex([c1,c2,c3], card);
    valList(idx) = tripletList(i).factorVal;
end

for i=1:n-2
    factors(i).var = [i, i+1, i+2];
    factors(i).card = card;
    factors(i).val = valList;
end