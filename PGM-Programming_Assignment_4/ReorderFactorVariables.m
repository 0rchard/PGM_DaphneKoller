function out = ReorderFactorVariables(in)

[S, I] = sort(in.var);

out.var = S;
out.card = in.card(I);

allAssignmentIn = IndexToAssignment(1:prod(in.card), in.card);
allAssignmentOut = allAssignmentIn(:,I);
out.val(AssignmentToIndex(allAssignmentOut, out.card)) = in.val;