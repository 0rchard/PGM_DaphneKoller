factors = repmat(struct('idx', [], 'val', []), 4, 1);
factors(1).idx = 1;
factors(1).val = 327.2;
factors(2).idx = 2;
factors(2).val = 368.2;
factors(3).idx = 3;
factors(3).val = 197.6;
factors(4).idx = 4;
factors(4).val = 178.4;

[unused, order] = sort([factors(:).val], 'descend');
sorted_fact = factors(order);