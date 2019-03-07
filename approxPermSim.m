function nu = approxPermSim(A, B)
  n = size(A, 1);

  U = kron(ones(n, 1), eye(n));
  V = kron(eye(n), ones(n, 1));

  u1 = U(:, 1);

  C = [U.' - ones(n, 1) * u1.'; V.' - ones(n, 1) * u1.'];

  C = [C(2:end, :); kron(B.', eye(n)) - kron(eye(n), A)]; % C is (2n + n^2) x n^2

  W = null(C); % TODO: what if C has full (column) rank?  Then W is empty...
  k = size(W, 2);

  [Q, R] = qr(W);

  E = [zeros(n^2 - k, k), eye(n^2 - k)];

  cvx_begin

  variable Z(n^2, n^2) semidefinite nonnegative

  maximize(trace(Z))

  subject to
  E * Q.' * Z * W == zeros(n^2 - k, k)    
  ones(1, n^2) * Z * ones(n^2, 1) == n^2
  
  cvx_end
  
  nu = cvx_optval;
end
