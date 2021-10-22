function sparse_amount = determineSparsityOfCTM(P)
% -------------------------------- %
% ---- Sparsity of CTM  ---------- %
% ---- @created: 2014-03-10 EK --- %
% ---- @revised: 2021-09-03 ND --- %
% -------------------------------- %

nCluster = size(P,1);
idx_IsZero = (P == 0);
sparse_amount = sum(sum(idx_IsZero))./(nCluster^2);

end
