function opt = sdpSpecialItr(A, B, i, j, minormax)
  n = size(A, 1);
  curr = zeros(n, n, 2);

  U = kron(ones(n, 1), eye(n));
  V = kron(eye(n), ones(n, 1));

  u1 = U(:, 1);

  C = [U.' - ones(n, 1) * u1.'; V.' - ones(n, 1) * u1.'];

  C = [C(2:end, :); kron(B.', eye(n)) - kron(eye(n), A)]; % C is (2n - 1 + n^2) x n^2
  
  % Orthonormalize the rows of C
  C = rref(C);
  k = rank(C);
  C = C(1 : k, :);

  tic;
  
  cvx_begin

  variable Z(n^2, n^2) nonnegative

  Xhat = diag(Z)
  X = reshape(Xhat, n, n)
  W = [Z, Xhat; Xhat', 1]

  if strcmp(minormax, 'min') || strcmp(minormax, 'minimum') || strcmp(minormax, 'minimize')
    minimize(X(i, j))
  else
    maximize(X(i, j))
  end

  subject to

  C * Z == 0;
  sum(sum(Z)) == n^2;
  
  X * ones(n, 1) == 1
  ones(1, n) * X == 1
  
  for u = 1 : n
    for v = 1 : n
      for w = 1 : n
	if u ~= w
	  Z(getIndex(u, v, n), getIndex(w, v, n)) == 0
	  Z(getIndex(v, u, n), getIndex(v, w, n)) == 0
	end
      end
    end
  end
  
  A * X == X * B

  for u = 1 : n
    for p = 1 : n
      sum(sum(Z(u : n : n^2, p : n : n^2) .* B)) == A(u, p)
      sum(sum(Z((u - 1) * n + 1 : u * n, (p - 1) * n + 1 : p * n) .* A)) == B(u, p)
    end
  end
  
  W == semidefinite(n^2 + 1)
  
  cvx_end

  toc;	
  disp(cvx_optval);

  opt = cvx_optval;

end
