function C = cycle(n)
  C = diag(ones(n - 1, 1), 1);
  C = C + C.';
  C(1, n) = 1;
  C(n, 1) = 1;
end
