function [L, U] = lpIter(A, B, L, U)
  n = size(A, 1);
  curr = zeros(n, n, 2);

  cvx_precision default
  
  for i = 1 : n
    for j = 1 : n
      for matr = 1 : 2
	cvx_begin quiet

	variable X(n, n) nonnegative

	if matr == 1
	  minimize(X(i, j))
	else
	  maximize(X(i, j))
	end

	subject to

	A * X == X * B
	X * ones(n, 1) == 1
	ones(1, n) * X == 1

	for p = 1 : n
	  for q = 1 : n
	    
	    if L(p, q) > 0 && U(p, q) < 1
	      U = Inf * ones(n, n);
	      L = -U;
	      return;
	    elseif L(p, q) > 0
	      X(p, q) == 1
	    elseif U(p, q) < 1
	      X(p, q) == 0
	    end
	    
	  end
	end
	
	cvx_end

	curr(i, j, matr) = cvx_optval;

      end
      
    end
  end

  L = curr(:, :, 1);
  U = curr(:, :, 2);
  
end
