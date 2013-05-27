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