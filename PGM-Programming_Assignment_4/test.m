t = struct('var', [], 'card', [], 'val', []);
t.var = C.nodes{1};
t.card = cardList(t.var);
[vars indexA] = unique([C.factorList.var], 'first');
idx_arr = zeros(1, length(vars));
K = length(C.factorList);
for i=1:length(vars)
    for k=1:K
        if idx_arr(i) == 0 && ~isempty(find(C.factorList(k).var == vars(i)));
            idx_arr(i) = k;
            break;
        end
    end
end

for v=1:length(t.var)
    idx = idx_arr(find(vars == t.var(v)));
    if v==1;
        factor = C.factorList(idx)
    else
        factor = FactorProduct(factor, C.factorList(idx))
    end
end