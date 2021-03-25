function Q = make_Q(n,eps,dist)
%MAKE_Q Summary of this function goes here
%   Detailed explanation goes here
    arguments
        n (1,1) int
        eps (1,1) double
        dist (1,1) string
    end

    switch dist
        case "grid"
            Q = make_Q_grid(n,eps);
        case "rand"
            err("Not implemented yet")
        otherwise
            err("dist not defined")
    end
end

function Q = make_Q_grid(n,eps)
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
        Q(i,i) = d*eps;
    end
end