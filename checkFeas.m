function feas = checkFeas(A, B, i, j, value)
  feas = true;

  n = size(A, 1);

  cvx_begin

  variable X(n, n) nonnegative

  maximize 1

  subject to

  A * X == X * B
  X * ones(n, 1) == 1
  ones(1, n) * X == 1
  X(i, j) == value
  
  cvx_end

  if ~strcmp(cvx_status, 'Solved')
    disp(cvx_status);
    feas = false;
  end
  
end
