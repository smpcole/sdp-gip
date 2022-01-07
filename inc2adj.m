function A = inc2adj(M)
  [n, m] = size(M);
  A = zeros(n, n);
  for j = 1 : m
    e = find(M(:, j));
    if length(e) == 1
      A(e(1), e(1)) = 1;
    elseif length(e) == 2
      A(e(1), e(2)) = 1;
      A(e(2), e(1)) = 1;
    end
  end
end
