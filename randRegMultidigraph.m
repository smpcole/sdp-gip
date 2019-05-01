function A = randRegMultidigraph(n, r)
  A = zeros(n, n);

  for i = 1 : r
    P = eye(n);
    A = A + P(randperm(n), :);
  end
  
end
