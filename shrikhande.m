function A = shrikhande(n)
    A = zeros(n^2, n^2);
    %T = zeros(n, n, n, n);
    for i = 1 : n
        for j = 1 : n
            k = i + 1;
            l = j + 1;
            if i == n
                k = 1;
            end
            if j == n
                l = 1;
            end
            topleft = (i - 1) * n + j;
            topright = (i - 1) * n + l;
            bottomleft = (k - 1) * n + j;
            bottomright = (k - 1) * n + l;
            vertices = [topleft, topright, bottomleft, bottomright];
            A(vertices, vertices) = 1 - blkdiag(1, ones(2, 2), 1);
%             A(topleft, topright) = 1;
%             A(topright, topleft) = 1;
%             A(topleft, bottomleft) = 1;
%             A(bottomleft, topleft) = 1;
%             A(topleft, bottomright) = 1;
%             A(bottomright, topleft) = 1;
%             A(topright, bottomright) = 1;
%             A(bottomright, topright) = 1;
%             A(bottomleft, bottomright) = 1;
%             A(bottomright, bottomleft) = 1;
%             T(i, j, i, l) = 1;
%             T(i, j, k, j) = 1;
%             T(i, j, k, l) = 1;
%             T(i, l, k, l) = 1;
%             T(k, j, k, l) = 1;
        end
    end
    
    
end

