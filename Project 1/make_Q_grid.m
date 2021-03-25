function Q = make_Q_grid(n,eps)
    arguments
        n (1,1) double {mustBeInteger}
        eps (1,1) double = 1/4
    end
    
    Q = zeros(n,n);
    for i = 1:25
        d = 0;
        if mod(i-1, 5) > 0
            Q(i, i-1) = eps;
            d = d+1;
        end
        if mod(i, 5) > 0
            Q(i, i+1) = eps;
            d = d+1;
        end
        if fix((i-1)/5) > 0
            Q(i, i-5) = eps;
            d = d+1;
        end
        if fix((i-1+5)/5) < 5
            Q(i, i+5) = eps;
            d = d+1;
        end
        Q(i,i) = 1-d*eps;
    end
end