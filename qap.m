function [y, Y] = qap(A, B, X0)

  uniqueopt = false;

  c = ceil(4 * (norm(A, 2) + norm(B, 2)));
  n = size(A, 1);
  Ac = A + c * eye(n);
  Bc = B + c * eye(n);

  fprintf('Initial pt:\n');
  disp(X0);

  while ~uniqueopt

    cvx_begin;

    variable Y(n, n) nonnegative;

    maximize trace(Ac * Y * Bc * X0' + Ac * X0 * Bc * Y');

    subject to;

    Y * ones(n, 1) == 1;
    ones(1, n) * Y == 1;
    
    cvx_end;

    fprintf('Initial LP:\n');
    disp(cvx_optval);
    disp(Y);

    y = cvx_optval;

    if all(rnd(Y) == rnd(X0))

      uniqueopt = true;

      fprintf('Checking if optimum is unique:\n');
      
      % Check if optimum is unique
      cvx_begin;

      variable Y(n, n) nonnegative;

      maximize(sum(sum(Y .* (X0 < .5))));

      subject to;

      Y * ones(n, 1) == 1;
      ones(1, n) * Y == 1;
      trace(Ac * Y * Bc * X0' + Ac * X0 * Bc * Y') >= y - .5;

      cvx_end;

      if cvx_optval > .5
	uniqueopt = false;
	fprintf('Found distinct optimum:\n');
	disp(Y);
	X0 = rnd(Y);
      end

    else
      X0 = rnd(Y);
    end
    
  end

  disp(trace(A * X0 * B * X0'));
  disp(sum(sum(A)));
  
end

function Y = rnd(X)
  Y = max(X > .75, X);
  Y = min(Y > .25, Y);
end
