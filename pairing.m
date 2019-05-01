% Generate a random k-regular graph on n vertices using the pairing method
function A = pairing(n, k)
  A = zeros(n, n);

  % Randomly permute the points; pairs of consecutive points will be matched
  pairs = randperm(n * k);

  % Keep track of each point's bucket.  There are n buckets with k points each
  buckets = kron(1 : n, ones(1, k));

  for p = 1 : 2 : (n * k)
    i = pairs(p);
    j = pairs(p + 1);

    % Add an edge from i's bucket to j's bucket
    u = buckets(i);
    v = buckets(j);
    A(u, v) = A(u, v) + 1;
    A(v, u) = A(u, v);
  end

  A = A + diag(diag(A)); % Loops contribute 2 toward the degree

end
