function [L, U, k] = lp(A, B)
  n = size(A, 1);
  L = zeros(n);
  U = ones(n);

  Lprev = L;
  Uprev = U;

  TOL = 10^-9;

  k = 0;
  
  while true
    [L, U] = lpIter(A, B, Lprev, Uprev);

    k = k + 1;
    
    if L(1, 1) == -Inf || U(1, 1) == Inf || all(all(abs([L - Lprev, U - Uprev]) < TOL))
      fprintf('Terminated after %d iteration(s)', k);
      break;
    end

    Lprev = L;
    Uprev = U;
    
  end
  
end
