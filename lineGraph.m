function L = lineGraph(G)
  n = size(G, 1);
  m = sum(sum(G)) / 2; % # of edges
  L = zeros(m, m);

  LtoG = zeros(m, 2);
  k = 0;

  for i = 1 : n
    for j = i : n
      for l = 1 : G(i, j)
	k = k + 1;
	LtoG(k, :) = [i, j];
      end
      
    end

  end

  disp(LtoG);

  assert(k == m);

  % For now assume no loops
  for i = 1 : m
    for j = i + 1 : m
      I = LtoG(i, :);
      J = LtoG(j, :);
      numEdges = length(intersect(I, J));
      L(i, j) = numEdges;
      L(j, i) = numEdges;
    end
  end
  
end
