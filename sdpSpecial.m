function [L, U] = sdpSpecial(A, B, savedata)
  n = size(A, 1);
  curr = zeros(n, n, 2);

  C = constraintMatrix(A, B);
  
  % Attempt to load saved data from file
  i = 1;
  j = 0;
  L = zeros(n);
  U = zeros(n);
  try
    load(savedata);
  catch err
    fprintf('Unable to load data from file %s\n', savedata);
  end
  
  j = j + 1; % Stored (i, j) are the last COMPLETED indices
  curr(:, :, 1) = L;
  curr(:, :, 2) = U;

  for i = i : n
    for j = j : n

      for matr = 1 : 2

	tic;
	
	cvx_begin

	variable Z(n^2, n^2) nonnegative

	Xhat = diag(Z)
	X = reshape(Xhat, n, n)
	W = [Z, Xhat; Xhat', 1]

	if matr == 1
	  minimize(X(i, j))
	else
	  maximize(X(i, j))
	end

	subject to

	C * Z == 0;
	sum(sum(Z)) == n^2;
	
	X * ones(n, 1) == 1
	ones(1, n) * X == 1
	
	for u = 1 : n
	  for v = 1 : n
	    for w = 1 : n
	      if u ~= w
		Z(getIndex(u, v, n), getIndex(w, v, n)) == 0
		Z(getIndex(v, u, n), getIndex(v, w, n)) == 0
	      end
	    end
	  end
	end
	
	A * X == X * B

	for u = 1 : n
	  for p = 1 : n
	    sum(sum(Z(u : n : n^2, p : n : n^2) .* B)) == A(u, p)
	    sum(sum(Z((u - 1) * n + 1 : u * n, (p - 1) * n + 1 : p * n) .* A)) == B(u, p)
	  end
	end
	
	W == semidefinite(n^2 + 1)
	
	cvx_end

	toc;	
	disp(cvx_optval);

	curr(i, j, matr) = cvx_optval;

      end
      
      L = curr(:, :, 1);
      U = curr(:, :, 2);

      try
	save(savedata, 'i', 'j', 'L', 'U');
      catch err
	fprintf('Warning: unable to save to file %s\n', savedata)
      end
      	
    end

    j = 1;

  end
  
end
