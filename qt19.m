function A = qt19
  K4 = ones(4, 4) - eye(4);
  A = kron(eye(3), K4);
  for i = 2 : 2 : 12
    j = mod(i - 1 + 3, 12) + 1;
    A(i, j) = 1;
    A(j, i) = 1;
  end
  
end
