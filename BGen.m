function B = BGen(A)
  n = size(A, 1);

  % Generate a random permutation matrix
  P = eye(n);
  P = P(randperm(n), :);

  B0 = P * A * P';

  X = randDblStoch(n, n^2, false);
  S = X - X';

  function Bt = gen(t)
    Bt = exp(t * S) * B0 * exp(-t * S);
  end

  B = @gen;
    
end
