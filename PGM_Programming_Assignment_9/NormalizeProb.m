function probN = NormalizeProb(prob)
denorm = sum(prob);
probN = prob ./ denorm;
