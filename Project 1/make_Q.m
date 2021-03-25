function Q = make_Q(n,dist)
%MAKE_Q Summary of this function goes here
%   Detailed explanation goes here
    arguments
        n (1,1) double {mustBeInteger}
        dist (1,1) string
    end

    switch dist
        case "grid"
            Q = make_Q_grid(n);
        case "rand"
            Q = make_Q_rand(n);
        otherwise
            err("dist not defined")
    end
end

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

function Q = make_Q_rand(n,r)
    arguments
        n (1,1) double {mustBeInteger}
        r (1,1) double = 4
    end
    
    x=10*rand(n, 1);
    y=10*rand(n, 1);
    
    max=0;
    for i=1:n
        din=0;
        for j=n:-1:i
            if norm([x(i), y(i)]-[x(j), y(j)])<=r && i~=j
                din=din+1;
            end
        end
        if din>max
            max=din;
        end
    end

    Q=eye(n);
    eps=1/max;

    for i=1:n
        for j=n:-1:i
            if norm([x(i), y(i)]-[x(j), y(j)])<=r && i~=j
                Q(i ,j)=eps;
                Q(i, i)=Q(i, i)-eps;
                Q(j, i)=eps;
                Q(j, j)=Q(j, j)-eps;
            end
        end
    end
end