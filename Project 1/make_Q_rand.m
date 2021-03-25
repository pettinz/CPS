function Q = make_Q_rand(n,x,y,r)
    arguments
        n (1,1) double {mustBeInteger}
        x (:,:) double
        y (:,:) double
        r (1,1) double = 4
    end
    
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