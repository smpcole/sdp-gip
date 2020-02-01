function nu = sdp4(A, B)
  n = size(A, 1);

  C = constraintMatrix(A, B);
  
  E = [kron(eye(n), ones(n, 1)), kron(ones(n, 1), eye(n))];
  
  cvx_begin

  variable Z(n^2, n^2) semidefinite nonnegative

  d = diag(Z)
  U = [Z, d; d', 1]
  
  maximize(trace(Z))

  subject to
  C * Z == zeros(k, n^2);

  Z * E == d * ones(1, 2 * n)
  d' * E == ones(1, 2 * n)

  for i = 1 : n
    for j = 1 : n
      sum(sum(Z(i : n : n^2, j : n : n^2) .* A)) == B(i, j)
      sum(sum(Z((i - 1) * n + 1 : i * n, (j - 1) * n + 1 : j * n) .* B)) == A(i, j)
    end
  end
  
  U == semidefinite(n^2 + 1)
  
  cvx_end
  
  nu = cvx_optval;
end
