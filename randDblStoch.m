function A = randDblStoch(n, k, symmetric)

  A = 0;

  w = exprnd(1, k, 1);
  w = w / sum(w);

  for i = 1 : k
    P = eye(n);
    P = P(randperm(n), :);
    A = A + w(i) * P;
  end

  if symmetric
    A = (A + A') / 2;
  end
  
  disp(A);

  if any(any(A <= 0))
    minpos = 1;
    for i = 1 : n
      for j = 1 : n
	if A(i, j) > 0 && A(i, j) < minpos
	  minpos = A(i, j);
	end
      end
    end
    t = minpos / 2;
    A = (1 - t) * A + t;
  end
    
end
