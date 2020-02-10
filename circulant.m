function A = circulant(n, l)
  A = zeros(n, n);
  for i = 1 : n
    for j = l
      A(i, mod(i - 1 + j, n) + 1) = 1;
      A(i, mod(i - 1 - j, n) + 1) = 1;
    end
  end
  
end
