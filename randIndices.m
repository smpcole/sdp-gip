function mask = randIndices(m, n, numindices, symm)

  max_indices = m * n;
  if symm
    assert(m == n);
    max_indices = n * (n + 1) / 2;
  end
  
  mask = false(m, n);
  
  if numindices > max_indices
    numindices = max_indices;
  end

  indices = randperm(max_indices);
  indices = sort(indices(1 : numindices));
  curr = 1;
  ij = 1;
  for i = 1 : m

    minj = 1;
    if symm
      minj = i;
    end

    for j = minj : n

      if curr <= numindices && ij == indices(curr)
	mask(i, j) = true;
	%if symm
	 % mask(j, i) = true;
	%end
	curr = curr + 1;
      end
      ij = ij + 1;
    end
  end
  
    
end
