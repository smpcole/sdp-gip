function [i, j] = vectorToMatrix(p, n)
  j = ceil(p / n);
  i = mod(p, n);
  if i == 0
    i = n;
  end
end
