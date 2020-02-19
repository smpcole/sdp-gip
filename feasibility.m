function feas = feasibility(A, B, num_nonneg, psd)
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


  model = ccp_model('sdp-gip');

  if psd
    Z = var_sdp(n^2, n^2);
  else
    Z = var_symm(n^2, n^2)
  end

  model.add_variable(Z);

  model.maximize(trace(Z));

  % C * Z == 0;
  for pq = 1 : n^2
    for i = 1 : r
      Cipq = zeros(n^2, n^2);
      Cipq(:, pq) = C(i, :);
      model.add_affine_constraint(inprod(Cipq, Z) == 0);
    end
  end

  model.add_affine_constraint(double(zeroindices) .* Z == zeros(n^2, n^2));

  % (2.7)
  for i = 1 : n
    ei = zeros(n, 1);
    ei(i) = 1;
    rowi = diag(kron(ones(n, 1), ei));
    coli = diag(kron(ei, ones(n, 1)));
    model.add_affine_constraint(inprod(rowi, Z) == 1);
    model.add_affine_constraint(inprod(coli, Z) == 1);
  end

  % Choose random indices to be nonnegative
  if num_nonneg >= n^2 * (n^2 - 1) / 2
    model.add_affine_constraint(Z >= 0);
  else

    nonneg_indices = double(randIndices(n^2, n^2, num_nonneg, true));
    model.add_affine_constraint(nonneg_indices .* Z >= 0);

  end
  
  model.solve;
  
end

function done
  fprintf('done\n');
end
