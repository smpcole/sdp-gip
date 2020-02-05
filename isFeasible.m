function feas = isFeasible(A, B, num_nonneg, psd, truncate)
  feas = true;

  if psd && truncate
    fprintf('Truncation is unsupported with positive-semidefinite constraint, so truncation will be ignored...\n');
    truncate = false;
  end

  n = size(A, 1);
  
  C = constraintMatrix(A, B);
  
  I = eye(n);
  J = ones(n, n);

  global zeroindices;
  zeroindices = kron(I, J - I) + kron(J - I, I); % (2.4)
  zeroindices = or(zeroindices, kron(J, A) ~= kron(B, J)); % (3.6)
  zeroindices = or(zeroindices, zeroindices');

  if truncate

    global Ztow;
    Ztow = zeros(n^2, n^2);
    global wtoZ;
    wtoZ = zeros(n^4, 2);
    curr = 1;
    for ij = 1 : n^2
      for pq = ij : n^2
	if ~zeroindices(ij, pq)
	  Ztow(ij, pq) = curr;
	  Ztow(pq, ij) = curr;
	  wtoZ(curr, :) = [ij, pq];
	  curr = curr + 1;
	end
      end
    end

    N = curr - 1;
    wtoZ = wtoZ(1 : N, :);

    global wcols;
    wcols = {};
    global Ccols;
    Ccols = {};
    for pq = 1 : n^2
      colmask = false(n^2, n^2);
      colmask(:, pq) = true;
      colmask = and(colmask, ~zeroindices);
      wcols{pq} = Ztow(colmask);
      Ccols{pq} = 1 : n^2;
      Ccols{pq} = Ccols{pq}(colmask(:, pq));
    end

    global wdiag;
    diagmask = and(eye(n^2), ~zeroindices);
    wdiag = Ztow(diagmask);

    wdiagrows = {};
    wdiagcols = {};
    for i = 1 : n
      rowmask = and(diag(mod(1 : n^2, n) == mod(i, n)), ~zeroindices);
      wdiagrows{i} = Ztow(rowmask);
    end
    for j = 1 : n
      colmask = false(1, n^2);
      colmask((j - 1) * n + 1 : j * n) = true;
      colmask = and(diag(colmask), ~zeroindices);
      wdiagcols{j} = Ztow(colmask);
    end

    cvx_begin

    variable w(N);
    display(N);
    display(length(w));

    subject to

    for pq = 1 : n^2
      C(:, Ccols{pq}) * w(wcols{pq}) == 0;
    end

    for i = 1 : n
      sum(w(wdiagrows{i})) == 1;
      sum(w(wdiagcols{i})) == 1;
    end

    if num_nonneg >= length(w)
      w >= 0;
    else
      % Choose random indices to be nonnegative
      nonneg_indices = randperm(N);
      nonneg_indices = nonneg_indices(1 : num_nonneg);
      w(nonneg_indices) >= 0;
    end
    
    cvx_end


  else

    cvx_begin

    variable Z(n^2, n^2) symmetric;

    Xhat = diag(Z);
    X = reshape(Xhat, n, n);
    W = [Z, Xhat; Xhat', 1];

    subject to

    C * Z == 0;

    Z .* zeroindices == 0; % (2.4) and (3.6)

    % (2.7)
    X * ones(n, 1) == 1;
    ones(1, n) * X == 1;

    if psd
      W == semidefinite(n^2 + 1);
    end

    % Choose random indices to be nonnegative

    if num_nonneg >= n^2 * (n^2 - 1) / 2
      Z >= 0;
    else

      nonneg_indices = randIndices(n^2, n^2, num_nonneg, true);
      Z(nonneg_indices) >= 0;

    end
    
    cvx_end;
    
  end  
  
  if cvx_optval == Inf
    feas = false;
  end
  
  
  
end

