function Q = make_Q_grid_uniform(n,eps)
    arguments
        n (1,1) double {mustBeInteger}
        eps (1,1) double = 0;
    end
    
    if eps == 0
        P = build_connections(n);
        eps = 1/(1+max(P(:)));
    end
    
    Q = zeros(n);
    for i = 1:n
        [nn, sn] = get_neighbors(i, n);
        for j = 1:sn
            Q(i,nn(j)) = eps;
        end
        Q(i,i) = 1-sn*eps;
    end
end