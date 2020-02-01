function nu = sdp3(A, B)
  n = size(A, 1);

  C = constraintMatrix(A, B);
  
  cvx_begin

  variable Z(n^2, n^2) semidefinite nonnegative

  d = diag(Z);
  U = [Z, d; d', 1];

  maximize(trace(Z))

  subject to
  U == semidefinite(n^2 + 1)
  C * Z == zeros(k, n^2);
  ones(1, n^2) * Z * ones(n^2, 1) == n^2
  
  cvx_end
  
  nu = cvx_optval;
end
