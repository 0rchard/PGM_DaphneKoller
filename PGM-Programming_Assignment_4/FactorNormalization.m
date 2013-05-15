function G = FactorNormalization(F)
G = F;
if isempty(G.val)
    return;
end
val_total = 0;
for i=1:prod(G.card)
    val_total = val_total + G.val(i);
end

G.val = G.val / val_total;
