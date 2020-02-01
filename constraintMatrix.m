function C = constraintMatrix(A, B)
  
  n = size(A, 1);
  
  U = kron(ones(n, 1), eye(n));
  V = kron(eye(n), ones(n, 1));

  u1 = U(:, 1);

  global C
  C = [U.' - ones(n, 1) * u1.'; V.' - ones(n, 1) * u1.'];

  C = [C(2:end, :); kron(B.', eye(n)) - kron(eye(n), A)]; % C is (2n - 1 + n^2) x n^2
  
  % Row reduce C
  C = rref(C);
  k = rank(C);
  C = C(1 : k, :);

end
