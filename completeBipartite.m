function A = completeBipartite(m, n)
    A = zeros(m + n, m + n);
    A(1 : m, m + 1 : end) = 1;
    A(m + 1 : end, 1 : m) = 1;
end

