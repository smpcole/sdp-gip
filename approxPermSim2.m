function nu = approxPermSim(A, B)
  n = size(A, 1);

  U = kron(ones(n, 1), eye(n));
  V = kron(eye(n), ones(n, 1));

  u1 = U(:, 1);

  C = [U.' - ones(n, 1) * u1.'; V.' - ones(n, 1) * u1.'];

  C = [C(2:end, :); kron(B.', eye(n)) - kron(eye(n), A)]; % C is (2n - 1 + n^2) x n^2

  % Orthonormalize the rows of C
  C = rref(C);
  k = rank(C);
  C = C(1 : k, :);
  
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
