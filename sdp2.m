function nu = sdp2(A, B)
  n = size(A, 1);

  C = constraintMatrix(A, B);
  
  cvx_begin

  variable U(n^2 + 1, n^2 + 1) semidefinite

  t = U(n^2 + 1, n^2 + 1)
  W = U(1 : n^2, 1 : n^2)
  Z = W - (t - 1) * eye(n^2)
  w = U(1 : n^2, n^2 + 1)
  d = diag(Z)

  maximize(trace(U) - (n^2 + 2) * t)

  subject to

  t >= 1
  t <= n + 2

  Z == semidefinite(n^2)
  Z >= 0
  C * Z == 0
  ones(1, n^2) * Z * ones(n^2, 1) == n^2

  w == d
  
  cvx_end

  fprintf('If the optimal value is < %d, then A and B are NOT permutation similar.\n', n - 1 - n^2);
  
  nu = cvx_optval;
end
