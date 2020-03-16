function feas = feasibility(A, B, num_nonneg)
  feas = true;

  n = size(A, 1);
  
  C = constraintMatrix(A, B);
  r = size(C, 1);
  
  I = eye(n);
  J = ones(n, n);

  global zeroindices;
  zeroindices = kron(I, J - I) + kron(J - I, I); % (2.4)
  zeroindices = or(zeroindices, kron(J, A) ~= kron(B, J)); % (3.6)
  zeroindices = or(zeroindices, zeroindices');
  zeroindices = triu(zeroindices); % Only care about upper triangular portion since it is symmetric

  blk{1, 1} = 's';
  blk{1, 2} = n^2;

  N = (n^2 + 2) * (n^2 + 1) / 2;
  
  % C * Z == 0;

  fprintf('Building null space constraints...');
  At{1} = sparse(N, r * n^2);
  next = 1;
  for i = 1 : r
    for pq = 1 : n^2
      Cipq = zeros(n^2 + 1, n^2 + 1);
      Cipq(1 : n^2, pq) = C(i, :) / 2;
      Cipq(pq, :) = Cipq(pq, :) + Cipq(:, pq)';
      At{1}(:, next) = svec(blk(1, :), Cipq);
      next = next + 1;
    end
  end
  done;

  b = zeros(r * n^2, 1);

  fprintf('Building zero index constraints...');
  numzeroindices = sum(sum(zeroindices));
  zeroconstraints = sparse(N, numzeroindices);
  next = 1;
  for ij = 1 : n^2
    for pq = ij : n^2
      if zeroindices(ij, pq)
	Eijpq = zeros(n^2 + 1, n^2 + 1);
	Eijpq(ij, pq) = .5;
	Eijpq(pq, ij) = Eijpq(pq, ij) + .5;
	zeroconstraints(:, next) = svec(blk(1, :), Eijpq);
	next = next + 1;
      end
    end
  end
  done;

  At{1} = [At{1}, zeroconstraints];
  b = [b; zeros(numzeroindices, 1)];

  fprintf('Building row/column sum constraints...');
  rowsums = sparse(N, n);
  colsums = sparse(N, n);;
  for i = 1 : n
    ei = zeros(n, 1);
    ei(i) = 1;
    rowi = diag([kron(ones(n, 1), ei); 0]);
    coli = diag([kron(ei, ones(n, 1)); 0]);
    rowsums(:, i) = svec(blk(1, :), rowi);
    colsums(:, i) = svec(blk(1, :), coli);
  end
  done;

  At{1} = [At{1}, rowsums, colsums];
  b = [b; ones(2 * n, 1)];

  fprintf('Building rank-1 PSD constraints...');
  Wconstr = sparse(N, n^2 + 1);
  for ij = 1 : n^2 + 1
    Wconstrij = zeros(n^2 + 1, n^2 + 1);
    Wconstrij(ij, ij) = 1;
    if ij <= n^2
      Wconstrij(ij, n^2 + 1) = -.5;
      Wconstrij(n^2 + 1, ij) = -.5;
    end
    Wconstr(:, ij) = svec(blk(1, :), Wconstrij);
  end
  b = [b; zeros(n^2, 1); 1];
  done;

  L = [];
  l = [];
  Bt = [];

  fprintf('Building nonnegative constraints...');
  % Choose random indices to be nonnegative
  if num_nonneg >= n^2 * (n^2 + 1) / 2
    L = 0;
  else

    nonneg_indices = double(randIndices(n^2, n^2, num_nonneg, true));

    l = zeros(num_nonneg, 1);

    Bt{1} = sparse(N, num_nonneg);

    next = 1;
    for ij = 1 : n^2
      for pq = ij : n^2
	if nonneg_indices(ij, pq)
	  Eijpq = zeros(n^2 + 1, n^2 + 1);
	  Eijpq(ij, pq) = .5;
	  Eijpq(pq, ij) = Eijpq(pq, ij) + .5;
	  Bt{1}(:, next) = svec(blk(1, :), Eijpq);
	end
      end
    end

  end
  done;

  objfun{1} = sparse(n^2, n^2);

  [obj, X, s, y, S, Z, ybar, v, info, runhist] = sdpnalplus(blk, At, objfun, b, L, [], Bt, l, []);
  
end

function done
  fprintf('done\n');
end
