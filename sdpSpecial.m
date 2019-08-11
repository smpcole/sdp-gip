function [L, U] = sdpSpecial(A, B)
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

  for i = 1 : n
    for j = 1 : n

      for matr = 1 : 2

	cvx_begin

	variable Z(n, n, n, n) nonnegative

	Zsq = reshape(Z, n^2, n^2)

	Xhat = diag(Zsq)
	X = reshape(Xhat, n, n)
	W = [Zsq, Xhat; Xhat', 1]

	if matr == 1
	  minimize(X(i, j))
	else
	  maximize(X(i, j))
	end

	subject to

	C * Zsq == 0;
	sum(sum(Zsq)) == n^2;
	
	X * ones(n, 1) == 1
	ones(1, n) * X == 1
	
	A * X == X * B

	for u = 1 : n
	  for p = 1 : n
	    sum(sum(reshape(Z(u, :, p, :), n, n) .* B)) == A(u, p)
	    sum(sum(reshape(Z(:, u, :, p), n, n) .* A)) == B(u, p)
	  end
	end
	
	W == semidefinite(n^2 + 1)
	
	cvx_end

	disp(cvx_optval);

	curr(i, j, matr) = cvx_optval;

      end

    end

  end

  L = curr(:, :, 1);
  U = curr(:, :, 2);
  
end