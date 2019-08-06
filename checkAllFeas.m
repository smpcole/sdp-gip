function [feas0, feas1] = checkAllFeas(A, B)

  n = size(A, 1);

  feas = true(n, n, 2);

  cvx_quiet true

  for i = 1 : n
    for j = 1 : n
      for value = 0 : 1
	feas(i, j, value + 1) = checkFeas(A, B, i, j, value);
      end
    end
  end
  
  feas0 = feas(:, :, 1);
  feas1 = feas(:, :, 2);
  
end
