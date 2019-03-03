function A = randSymmStochastic(n)
  A = zeros(n, n);

  for i = 1 : n
    len = n - i + 1;
    rem = 1 - sum(A(i, :));
    row = rand(1, len);
    row = row / sum(row) * rem;
    A(i, i : end) = row;
    A(i : end, i) = row.';
  end
  
end
